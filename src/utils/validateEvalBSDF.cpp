#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/texture.h>
#include <boost/filesystem/path.hpp>

MTS_NAMESPACE_BEGIN

class ValidateEvalBSDF : public Utility {
public:
	typedef std::vector<std::vector<Vector3d> > Vector2DArray;
	typedef std::vector<std::vector<Point2> > Sample2DArray;

	int run(int argc, char **argv) {
		m_wi = normalize(Vector(std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3])));

		m_size = std::atoi(argv[4]);
		m_sqrtSpp = std::atoi(argv[5]);
		m_spp = m_sqrtSpp * m_sqrtSpp;

		char *outFilename = argv[6];

		//m_wi = normalize(Vector(0.3, 0.2, 1));
		//m_wi = normalize(Vector(0, 0, 1));

		Properties props;
		props = Properties("independent");
		props.setInteger("seed", 19931004);
		m_sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), props));
		m_sampler->configure();

		// generate samples
		std::vector<std::vector<Point2> > samples(m_sqrtSpp, std::vector<Point2>(m_sqrtSpp));
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				samples[i][j] = m_sampler->next2D();
			}
		}

		Vector2DArray bsdfValues(m_size, std::vector<Vector3d>(m_size));

		/*
		props = Properties("aniso_roughdiffuse_simple");
		props.setSpectrum("reflectance", Spectrum(1.0f));
		props.setFloat("alpha", 1.0);
		BSDF *bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props));
		bsdf->configure();

		calcEvalBSDF(bsdf, m_wi, samples, bsdfValues);
		outputBitmap(bsdfValues, "microfacet_eval.exr");

		props = Properties("aniso_roughdiffuse");
		props.setSpectrum("reflectance", Spectrum(1.0f));
		props.setBoolean("sampleVisibility", true);
		bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props));
		bsdf->configure();

		calcEvalBSDF(bsdf, m_wi, samples, bsdfValues);
		outputBitmap(bsdfValues, "microfacet_LEADR_eval.exr");
		*/

		/*
		props = Properties("roughdiffuse_multi");
		props.setSpectrum("reflectance", Spectrum(0.9f));
		props.setFloat("alphaU", 0.5);
		props.setFloat("alphaV", 0.5);
		props.setInteger("scatteringOrderMax", 2);
		BSDF *bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props));
		bsdf->configure();
		*/

		/*
		props = Properties("diffuse");
		props.setSpectrum("reflectance", Spectrum(0.9f));
		BSDF *baseBSDF = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props));
		baseBSDF->configure();

		props = Properties("roughbsdf_multi");
		props.setInteger("scatteringOrderMax", 10);
		Spectrum moments1;
		moments1[0] = 5.125; moments1[1] = 5.125; moments1[2] = 5.0;
		props.setSpectrum("moments1", moments1);
		BSDF *bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props));
		bsdf->addChild("baseBSDF", baseBSDF);
		bsdf->configure();
		*/

		props = Properties("roughconductor");
		props.setString("distribution", "GGX");
		props.setFloat("alpha", 0.1);
		props.setString("material", "none");
		Spectrum spec;
		spec[0] = 0.9; spec[1] = 0.1; spec[2] = 0.1;
		props.setSpectrum("specularReflectance", spec);
		
		BSDF *bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props));
		bsdf->configure();

		/*
		Intersection its;
		its.baseFrame = Frame(Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1));
		its.shFrame = its.baseFrame;
		BSDFSamplingRecord bRec(its, m_wi, Vector(0, 0, 1));
		Spectrum tmp = bsdf->eval(bRec);
		Log(EInfo, "eval res = (%.6f, %.6f, %.6f)", tmp[0], tmp[1], tmp[2]);
		*/

		calcEvalBSDF(bsdf, m_wi, samples, bsdfValues);
		outputBitmap(bsdfValues, outFilename);

		return 0;
	}

	void calcEvalBSDF(BSDF *bsdf, const Vector &wi, const Sample2DArray &samples, Vector2DArray &bsdfValues) {
		for (int i = 0; i < m_size; i++)
			for (int j = 0; j < m_size; j++)
				bsdfValues[i][j] = Vector3d(0.0f);

#pragma omp parallel for
		for (int i = 0; i < m_size; i++) {
			for (int j = 0; j < m_size; j++) {
				Vector3d res(0.0);
				for (int x = 0; x < m_sqrtSpp; x++) {
					for (int y = 0; y < m_sqrtSpp; y++) {
						double dx = (j + samples[x][y].x) / (double)m_size * 2.0 - 1.0;
						double dy = (i + samples[x][y].y) / (double)m_size * 2.0 - 1.0;
						if (dx * dx + dy * dy > 1.0)
							continue;
						double dz = sqrt(1.0 - dx * dx - dy * dy);

						Intersection its;
						its.baseFrame = Frame(Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1));
						its.shFrame = its.baseFrame;
						BSDFSamplingRecord bRec(its, wi, Vector(dx, dy, dz));
						Spectrum tmp = bsdf->eval(bRec);

						for (int c = 0; c < 3; c++) {
							res[c] += tmp[c];
						}
					}
				}

				res /= (double)m_spp;

				for (int c = 0; c < 3; c++) {
					bsdfValues[m_size - i - 1][j][c] = res[c];
				}
			}
		}
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

	ref<Sampler> m_sampler;
	int m_size;
	int m_sqrtSpp, m_spp;
	Vector3 m_wi;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(ValidateEvalBSDF, "Visualize BSDF eval")
MTS_NAMESPACE_END