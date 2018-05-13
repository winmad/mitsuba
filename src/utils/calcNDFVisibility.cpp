#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <boost/filesystem/path.hpp>
#include "../bsdfs/microfacet.h"
#include "multi_vmf.h"

MTS_NAMESPACE_BEGIN

class CalcNDFVisibility : public Utility {
public:
	int run(int argc, char **argv) {
		m_mode = std::atoi(argv[1]);
		char *inFilename = argv[2];
		m_size = std::atoi(argv[3]);
		char *outFilename = argv[4];
		
		m_res = new Bitmap(Bitmap::ELuminance, Bitmap::EFloat32, Vector2i(m_size, m_size));

		if (m_mode == 1) {
			fromTabulatedNDF(inFilename, outFilename);
		} else if (m_mode == 2) {
			fromVmfNDF(inFilename, outFilename);
		}

		return 0;
	}

	inline bool getDir(int r, int c, int size, Vector &d) {
		d.x = (c + 0.5) / (double)size * 2.0 - 1.0;
		d.y = ((size - r - 1) + 0.5) / (double)size * 2.0 - 1.0;
		if (d.x * d.x + d.y * d.y >= 1.0)
			return false;
		d.z = std::sqrt(1.0 - d.x * d.x - d.y * d.y);
		return true;
	}

	Float G1(const Vector &wi, const Vector &nMeso, const Vector2i &ndfSize, float *data) {
		Float integral = 0.0;
		for (int r2 = 0; r2 < ndfSize.y; r2++) {
			for (int c2 = 0; c2 < ndfSize.x; c2++) {
				Vector wm;
				if (!getDir(r2, c2, ndfSize.x, wm))
					continue;

				Float dwm = 1.0 / wm.z * 4.0 / (ndfSize.x * ndfSize.y);
				Float value = data[r2 * ndfSize.x + c2] / wm.z;
				integral += std::max(0.0, dot(wi, wm)) * value * dwm;
			}
		}

		Float res = std::max(0.0, dot(wi, nMeso)) / std::max(0.0, nMeso.z);

		//Log(EInfo, "\n=============\nwi = (%.6f, %.6f, %.6f)\nnMeso = (%.6f, %.6f, %.6f)\n%.6f / %.6f",
		//	wi.x, wi.y, wi.z, nMeso.x, nMeso.y, nMeso.z, res, integral);

		res /= integral;
		return res;
	}

	void fromTabulatedNDF(char *ndfFilename, char *outFilename) {
		ref<Bitmap> ndf = new Bitmap(fs::path(ndfFilename));
		Vector2i ndfSize = ndf->getSize();
		float *data = ndf->getFloat32Data();
		
		Vector nMeso(0.0);
		for (int r = 0; r < ndfSize.y; r++) {
			for (int c = 0; c < ndfSize.x; c++) {
				Float value = data[r * ndfSize.x + c];
				Vector d;
				if (getDir(r, c, ndfSize.x, d)) {
					nMeso += value * d;
				}
			}
		}
		Float len = nMeso.length();
		if (len > 1e-4)
			nMeso /= len;

		float *res = m_res->getFloat32Data();
#pragma omp parallel for
		for (int r = 0; r < m_size; r++) {
			for (int c = 0; c < m_size; c++) {
				Vector wi;
				if (!getDir(r, c, m_size, wi))
					continue;

				res[r * m_size + c] = G1(wi, nMeso, ndfSize, data);
			}
		}

		ref<FileStream> stream = new FileStream(outFilename, FileStream::ETruncWrite);
		m_res->write(Bitmap::EOpenEXR, stream);
	}

	inline Float conv(Float w1, Float lambda1, const Vector &mu1,
			Float w2, Float lambda2, const Vector &mu2) const {
		Vector d = mu1 * lambda1 + mu2 * lambda2;
		Float len = d.length();
		Float res = w1 * w2 * 4.0 * M_PI * math::fastexp(-lambda1 - lambda2) 
			* std::sinh(len) / len;
		return res;
	}

	Float G1(const Vector &wi, const Vector &nMeso, const std::vector<Float> &alpha,
			const std::vector<Float> &kappa, const std::vector<Vector> &mu) {
		Float wCos = 1.1767;
		Float lambdaCos = 2.1440;
		Float integral = 0.0;
		
		Float vertProjArea = 0.0;
		Vector wg(0, 0, 1);
		
		for (int i = 0; i < alpha.size(); i++) {
			if (alpha[i] < 1e-4)
				continue;
			Float wNDF = alpha[i] * kappa[i] / (2 * M_PI * (1 - math::fastexp(-2 * kappa[i])));
			integral += conv(wCos, lambdaCos, wi, wNDF, kappa[i], mu[i]);
			vertProjArea += conv(wCos, lambdaCos, wg, wNDF, kappa[i], mu[i]);
		}

		Float res = std::max(0.0, dot(wi, nMeso)) / std::max(0.0, nMeso.z);
		integral /= vertProjArea;

		//Log(EInfo, "\n=============\nwi = (%.6f, %.6f, %.6f)\nnMeso = (%.6f, %.6f, %.6f)\n%.6f / %.6f",
		//	wi.x, wi.y, wi.z, nMeso.x, nMeso.y, nMeso.z, res, integral);

		res /= integral;
		return res;
	}

	void fromVmfNDF(char *paramFilename, char *outFilename) {
		FILE *fp = fopen(paramFilename, "r");
		int numLobes;
		fscanf(fp, "%d", &numLobes);
		std::vector<Float> alpha(numLobes);
		std::vector<Float> kappa(numLobes);
		std::vector<Vector> mu(numLobes);
		for (int i = 0; i < numLobes; i++) {
			Float x, y, z;
			fscanf(fp, "%lf %lf %lf %lf %lf", &alpha[i], &kappa[i], &x, &y, &z);
			mu[i].x = x;
			mu[i].y = y;
			mu[i].z = z;
		}

		Vector nMeso(0.0);
		for (int i = 0; i < numLobes; i++) {
			nMeso += alpha[i] * mu[i];
		}
		Float len = nMeso.length();
		if (len > 1e-4)
			nMeso /= len;

		float *res = m_res->getFloat32Data();
#pragma omp parallel for
		for (int r = 0; r < m_size; r++) {
			for (int c = 0; c < m_size; c++) {
				Vector wi;
				if (!getDir(r, c, m_size, wi))
					continue;

				res[r * m_size + c] = G1(wi, nMeso, alpha, kappa, mu);
			}
		}

		ref<FileStream> stream = new FileStream(outFilename, FileStream::ETruncWrite);
		m_res->write(Bitmap::EOpenEXR, stream);
	}

	int m_mode;
	ref_vector<Sampler> m_samplers;
	int m_size;
	ref<Bitmap> m_res;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(CalcNDFVisibility, "Compute the spherical visibility function of a patch in the heightmap")
MTS_NAMESPACE_END