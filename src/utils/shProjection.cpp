#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <boost/filesystem/path.hpp>

MTS_NAMESPACE_BEGIN

struct BRDFFunctor {
	BRDFFunctor(BSDF *bsdf, const Vector &wi, int channel) {
		m_bsdf = bsdf;
		m_wi = wi;
		m_channel = channel;
	}

	Float operator()(const Vector &w) const {
		Intersection its;
		BSDFSamplingRecord bRec(its, m_wi, w);
		Spectrum res = m_bsdf->eval(bRec);
		return res[m_channel];
	}

	ref<BSDF> m_bsdf;
	Vector m_wi;
	int m_channel;
};

struct LobeFunctor {
	LobeFunctor(char *filename, const Vector &wi, int channel) {
		m_bitmap = new Bitmap(fs::path(filename));
		m_wi = wi;
		m_channel = channel;
	}

	Float operator()(int r, int c) const {
		return m_bitmap->getPixel(Point2i(c, r))[m_channel];
	}

	int getRes() const {
		return m_bitmap->getSize()[0];
	}

	ref<Bitmap> m_bitmap;
	Vector m_wi;
	int m_channel;
};

class SHProjection : public Utility {
public:
	int run(int argc, char **argv) {
		m_bands = std::atoi(argv[1]);

		/*
		m_resolution = std::atoi(argv[2]);
		m_alpha = std::atof(argv[3]);

		Properties props = Properties("roughconductor");
		props.setString("distribution", "beckmann");
		props.setFloat("alpha", m_alpha);
		props.setString("material", "Al");
		BSDF *bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props));
		bsdf->configure();

		BRDFFunctor f(bsdf, Vector(0, 0, 1), 0);
		SHVector sh(m_bands);
		sh.project(f, m_resolution);
		
		Log(EInfo, "bands = %d", m_bands);
		for (int i = 0; i < m_bands; i++) {
			Log(EInfo, "%.8f", sh(i, 0));
		}

		outputBRDFLobe(bsdf);
		*/

		LobeFunctor fLobe(argv[2], Vector(0, 0, 1), 0);
		SHVector shLobe(m_bands);
		shLobe.projectLobe(fLobe, fLobe.getRes());
		
		Log(EInfo, "bands = %d", m_bands);
		for (int i = 0; i < m_bands; i++) {
			std::ostringstream oss;
			for (int j = -i; j <= i; j++) {
				oss << shLobe(i, j) << " ";
			}
			Log(EInfo, "%s", oss.str().c_str());
		}
		
		Log(EInfo, "====== Energy =======");
		std::ostringstream oss;
		for (int i = 0; i < m_bands; i++)
			oss << shLobe.energy(i) << " ";
		Log(EInfo, "%s", oss.str().c_str());

		outputSHevalLobe(shLobe, argv[3]);

		return 0;
	}

	/// Helper function: creates 1D uniform random numbers
	inline Float rand1D() const {
		return std::min(0.9999, Float(rand()) / Float(RAND_MAX));
	}

	/// Helper function: creates 2D uniform random numbers
	inline Point2 rand2D() const {
		return Point2(rand1D(), rand1D());
	}

	void outputBRDFLobe(BSDF *bsdf) {
		int m_size = 128;
		std::vector<std::vector<Vector3d> > bsdfValues(m_size, std::vector<Vector3d>(m_size));
		for (int i = 0; i < m_size; i++)
			for (int j = 0; j < m_size; j++)
				bsdfValues[i][j] = Vector3d(0.0f);

		int sampleCount = 100;
#pragma omp parallel for
		for (int i = 0; i < m_size; i++) {
			for (int j = 0; j < m_size; j++) {
				Vector3d res(0.0);
				for (int k = 0; k < sampleCount; k++) {
					double x = (j + rand1D()) / (double)m_size * 2.0 - 1.0;
					double y = (i + rand1D()) / (double)m_size * 2.0 - 1.0;
					if (x * x + y * y > 1.0)
						continue;
					double z = sqrt(1.0 - x * x - y * y);

					Intersection its;
					BSDFSamplingRecord bRec(its, Vector(0, 0, 1), Vector(x, y, z));
					Spectrum tmp = bsdf->eval(bRec);

					for (int c = 0; c < 3; c++) {
						res[c] += tmp[c];
					}
				}

				res /= (double)sampleCount;

				for (int c = 0; c < 3; c++) {
					bsdfValues[m_size - i - 1][j][c] = res[c];
				}
			}
		}

		ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, Vector2i(m_size));
		float *data = bitmap->getFloat32Data();
		for (int r = 0; r < m_size; r++) {
			for (int c = 0; c < m_size; c++) {
				*data++ = bsdfValues[r][c][0];
				*data++ = bsdfValues[r][c][1];
				*data++ = bsdfValues[r][c][2];
			}
		}
		fs::path filename("lobe_brdf.exr");
		filename.replace_extension(".exr");
		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		bitmap->write(Bitmap::EOpenEXR, stream);
	}

	void outputSHevalLobe(SHVector &sh, char *filename) {
		int m_size = 128;
		std::vector<std::vector<Vector3d> > bsdfValues(m_size, std::vector<Vector3d>(m_size));
		for (int i = 0; i < m_size; i++)
			for (int j = 0; j < m_size; j++)
				bsdfValues[i][j] = Vector3d(0.0f);

		int sampleCount = 100;
#pragma omp parallel for
		for (int i = 0; i < m_size; i++) {
			for (int j = 0; j < m_size; j++) {
				Vector3d res(0.0);
				for (int k = 0; k < sampleCount; k++) {
					double x = (j + rand1D()) / (double)m_size * 2.0 - 1.0;
					double y = (i + rand1D()) / (double)m_size * 2.0 - 1.0;
					if (x * x + y * y > 1.0)
						continue;
					double z = sqrt(1.0 - x * x - y * y);

					//Spectrum tmp = bsdf->eval(bRec);
					Float tmp = sh.eval(Vector(x, y, z));

					for (int c = 0; c < 3; c++) {
						res[c] += tmp;
					}
				}

				res /= (double)sampleCount;

				for (int c = 0; c < 3; c++) {
					bsdfValues[m_size - i - 1][j][c] = res[c];
				}
			}
		}

		ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, Vector2i(m_size));
		float *data = bitmap->getFloat32Data();
		for (int r = 0; r < m_size; r++) {
			for (int c = 0; c < m_size; c++) {
				*data++ = bsdfValues[r][c][0];
				*data++ = bsdfValues[r][c][1];
				*data++ = bsdfValues[r][c][2];
			}
		}

		ref<FileStream> stream = new FileStream(fs::path(filename), FileStream::ETruncWrite);
		bitmap->write(Bitmap::EOpenEXR, stream);
	}

	int m_bands;
	int m_resolution;
	Float m_alpha;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(SHProjection, "Project spherical function to SH basis")
MTS_NAMESPACE_END