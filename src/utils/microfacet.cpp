#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/core/vmf.h>
#include <mitsuba/core/warp.h>
#include "../bsdfs/microfacet.h"

MTS_NAMESPACE_BEGIN

class TestMicrofacet : public Utility {
public:
	typedef std::vector<std::vector<Vector3d> > Vector2DArray;

	int run(int argc, char **argv) {
		MicrofacetDistribution dist(MicrofacetDistribution::EBeckmann, std::atof(argv[1]));
		m_size = std::atoi(argv[2]);
		char *filename = argv[3];

		Vector2DArray data;
		data.resize(m_size);
		for (int i = 0; i < m_size; i++)
			data[i].resize(m_size);

		for (int r = 0; r < m_size; r++) {
			for (int c = 0; c < m_size; c++) {
				Point2 p((c + 0.5) / m_size, (r + 0.5) / m_size);
				Vector w = warp::squareToUniformHemisphereConcentric(p);
				Float val = dist.eval(w);
				for (int k = 0; k < 3; k++)
					data[r][c][k] = val;
			}
		}

		double dp = 2 * M_PI / (double)(m_size * m_size);
		double totProjArea = 0.0;
		for (int r = 0; r < m_size; r++) {
			double y = (r + 0.5) / (double)m_size;
			for (int c = 0; c < m_size; c++) {
				double x = (c + 0.5) / (double)m_size;
				Vector normal = warp::squareToUniformHemisphereConcentric(Point2(x, y));
				totProjArea += data[r][c][0] * normal.z * dp;
			}
		}
		Log(EInfo, "Validation: int[D(wm)cos(wm)] = %.8f should be equal to 1", totProjArea);

		outputBitmap(data, filename);
		
		/*		
		Vector m(0, 0, 1);
		
		Vector wi(std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]));
		wi = normalize(wi);

		Log(EInfo, "G1 = %.6f", dist.smithG1(wi, m));

		int nTheta = 180;
		int nPhi = 360;
		double dTheta = 0.5 * M_PI / (double)nTheta;
		double dPhi = 2.0 * M_PI / (double)nPhi;

		double totVisible = 0.0;
		for (int i = 0; i < nTheta; i++) {
			double theta = (i + 0.5f) * dTheta;
			double sinTheta = std::sin(theta);
			for (int j = 0; j < nPhi; j++) {
				double phi = (j + 0.5f) * dPhi;
				Vector wm(sinTheta * std::cos(phi), sinTheta * std::sin(phi), std::cos(theta));
				double cosTerm = std::max(0.0, dot(wi, wm));
				totVisible += dist.eval(wm) * cosTerm * sinTheta * dTheta * dPhi;
			}
		}
		Log(EInfo, "total visible = %.6f", totVisible);

		double res = 0.0;
		for (int i = 0; i < nTheta; i++) {
			double theta = (i + 0.5f) * dTheta;
			double sinTheta = std::sin(theta);
			for (int j = 0; j < nPhi; j++) {
				double phi = (j + 0.5f) * dPhi;
				Vector wo(sinTheta * std::cos(phi), sinTheta * std::sin(phi), std::cos(theta));
				Vector wh = wi + wo;
				wh = normalize(wh);
				double cosTerm = std::max(0.0, dot(wi, wh));
				double visD_wh = dist.eval(wh) * cosTerm / totVisible;
				//double shadowGivenMaskTerm = (1.0 + dist.smithLambda(wi, m)) / (1.0 + dist.smithLambda(wi, m) + dist.smithLambda(wo, m));
				double shadowGivenMaskTerm = dist.smithG1(wo, m);
				res += shadowGivenMaskTerm * visD_wh * sinTheta * dTheta * dPhi / (4.0 * wh.z);
			}
		}
		
		Log(EInfo, "avg G2 = %.6f %.6f", res, res * dist.smithG1(wi, m));
		*/
		return 0;
	}

	void outputBitmap(const Vector2DArray &bsdfValues, char *filename) const {
		ref<Bitmap> bitmap = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, Vector2i(m_size));
		float *data = bitmap->getFloat32Data();
		for (int r = 0; r < m_size; r++) {
			for (int c = 0; c < m_size; c++) {
				*data++ = bsdfValues[r][c][0];
				*data++ = bsdfValues[r][c][1];
				*data++ = bsdfValues[r][c][2];
			}
		}
		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		bitmap->write(Bitmap::EOpenEXR, stream);
	}

	int m_size;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(TestMicrofacet, "Play with microfacet distribution")
MTS_NAMESPACE_END
