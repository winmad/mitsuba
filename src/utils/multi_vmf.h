#if !defined(__MULTI_VMF_H)
#define __MULTI_VMF_H

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/vmf.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <boost/filesystem/path.hpp>

MTS_NAMESPACE_BEGIN

struct MultiLobeVMF {
	MultiLobeVMF(int numLobes = 0) : m_numLobes(numLobes) {
		m_alpha.resize(numLobes);
		m_mu.resize(numLobes);
		m_dist.resize(numLobes);
	}

	Float eval(const Vector &w) const {
		Float res = 0.0;
		for (int i = 0; i < m_numLobes; i++) {
			if (m_alpha[i] < 1e-8)
				continue;
			Float cosTheta = dot(w, m_mu[i]);
			res += m_dist[i].eval(cosTheta) * m_alpha[i];
		}
		return res;
	}

	/// Helper function: creates 1D uniform random numbers
	inline Float rand1D() const {
		return std::min(0.9999, Float(rand()) / Float(RAND_MAX));
	}

	void outputDistribution(int size, char *filename) {
		std::vector<std::vector<Vector3d> > values(size, std::vector<Vector3d>(size));
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				values[i][j] = Vector3d(0.0f);

		// ensure \int D_{wi}(wm)dwm = 1
		// integrate over the whole sphere
		double dp = 2 * M_PI / (double)(size * size);
		double totD = 0.0;

		int sampleCount = 1;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				Point2 p;
				p.x = (j + 0.5) / (double)size;
				p.y = (i + 0.5) / (double)size;
				Vector w = warp::squareToUniformHemisphereConcentric(p);
				Float tmp = eval(w);

				for (int c = 0; c < 3; c++) {
					values[i][j][c] = tmp;
				}

				totD += tmp * dp;

				// the other hemisphere
				w.z = -w.z;
				tmp = eval(w);
				totD += tmp * dp;
			}
		}

		ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, Vector2i(size));
		float *data = bitmap->getFloat32Data();
		for (int r = 0; r < size; r++) {
			for (int c = 0; c < size; c++) {
				*data++ = values[r][c][0];
				*data++ = values[r][c][1];
				*data++ = values[r][c][2];
			}
		}

		ref<FileStream> stream = new FileStream(fs::path(filename), FileStream::ETruncWrite);
		bitmap->write(Bitmap::EOpenEXR, stream);

		printf("Validation: int[D_wi(wm)] = %.8f should be equal to 1\n", totD);
	}

	int m_numLobes;
	std::vector<Float> m_alpha;
	std::vector<Vector> m_mu;
	std::vector<VonMisesFisherDistr> m_dist;
};

MTS_NAMESPACE_END

#endif