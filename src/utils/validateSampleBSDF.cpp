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

class ValidateSampleBSDF : public Utility {
public:
	typedef std::vector<std::vector<Vector3d> > Vector2DArray;
	typedef std::vector<std::vector<Point2> > Sample2DArray;

	int run(int argc, char **argv) {
		m_size = std::atoi(argv[1]);
		m_sqrtSpp = std::atoi(argv[2]);
		m_spp = m_sqrtSpp * m_sqrtSpp;

		m_wi = normalize(Vector(0.0, 0.0, 1.0));

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
		props = Properties("diffuse");
		props.setSpectrum("reflectance", Spectrum(1.0f));
		*/

		props = Properties("roughconductor");
		props.setString("distribution", "GGX");
		props.setFloat("alpha", 0.1);
		props.setString("material", "none");
		Spectrum spec;
		spec[0] = 0.9; spec[1] = 0.1; spec[2] = 0.1;
		props.setSpectrum("specularReflectance", spec);

		BSDF *baseBSDF = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props));
		baseBSDF->configure();

		props = Properties("multilobe_bsdf");
		props.setInteger("numLobes", 6);
		props.setString("prefix", "F:\\heightfield_prefiltering\\data\\downsample_gabardine\\glossy_env\\vmf_8x_lobe_");
		props.setFloat("uvscale", 4);
		BSDF *bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props));
		bsdf->addChild("baseBSDF", baseBSDF);
		bsdf->configure();

		/*
		props = Properties("sh_scaled_bsdf");
		props.setInteger("numCoeffs", 9);
		props.setString("prefix", "F:\\heightfield_prefiltering\\data\\downsample_gabardine\\glossy_env\\albedo_8x_direct_sh_");
		props.setFloat("uvscale", 4);
		BSDF *scaledBSDF = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props));
		scaledBSDF->addChild(bsdf);
		scaledBSDF->configure();
		*/

		props = Properties("tabulated_scaled_bsdf");
		props.setString("angularScaleFilename", 
			"F:\\heightfield_prefiltering\\data\\downsample_gabardine\\glossy_change_view\\result_rLobes_8x_blocksize_1024\\A_matrix.exr");
		BSDF *scaledBSDF = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props));
		scaledBSDF->addChild(bsdf);
		scaledBSDF->configure();

		Intersection its;
		its.shFrame = Frame(Vector(0, 0, 1));
		its.geoFrame = its.shFrame;
		its.baseFrame = its.shFrame;
		its.uv = Point2(0.5, 0.5);
// 		BSDFSamplingRecord bRec(its, Vector(0, 0, 1), Vector(std::sqrt(0.5), 0.0, std::sqrt(0.5)));
// 		Spectrum tmp = bsdf->eval(bRec);
		
// 		its.wi = Vector(0, 0, 1);
// 		BSDFSamplingRecord bRec(its, m_sampler);
// 		Spectrum tmp = bsdf->sample(bRec, m_sampler->next2D());
// 		Log(EInfo, "w = (%.6f, %.6f, %.6f), wo = (%.6f, %.6f, %.6f)",
// 			tmp[0], tmp[1], tmp[2], bRec.wo.x, bRec.wo.y, bRec.wo.z);

		calcEvalBSDF(bsdf, m_wi, samples, bsdfValues);
		outputBitmap(bsdfValues, "multilobe_bsdf_eval.exr");
 
		calcSampleBSDF(bsdf, m_wi, bsdfValues);
		outputBitmap(bsdfValues, "multilobe_bsdf_sample.exr");

		calcEvalBSDF(scaledBSDF, m_wi, samples, bsdfValues);
		outputBitmap(bsdfValues, "scaled_multilobe_bsdf_eval.exr");

		calcSampleBSDF(scaledBSDF, m_wi, bsdfValues);
		outputBitmap(bsdfValues, "scaled_multilobe_bsdf_sample.exr");
		
		return 0;
	}

	void calcEvalBSDF(BSDF *bsdf, const Vector &wi, const Sample2DArray &samples, Vector2DArray &bsdfValues) {
		for (int i = 0; i < m_size; i++)
			for (int j = 0; j < m_size; j++)
				bsdfValues[i][j] = Vector3d(0.0f);

		Intersection its;
		its.shFrame = Frame(Vector(0, 0, 1));
		its.geoFrame = its.shFrame;
		its.baseFrame = its.shFrame;
		its.uv = Point2(0.5, 0.5);

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

	void calcSampleBSDF(BSDF *bsdf, const Vector &wi, Vector2DArray &bsdfValues) {
		for (int i = 0; i < m_size; i++)
			for (int j = 0; j < m_size; j++)
				bsdfValues[i][j] = Vector3d(0.0f);

		Intersection its;
		its.shFrame = Frame(Vector(0, 0, 1));
		its.geoFrame = its.shFrame;
		its.baseFrame = its.shFrame;
		its.uv = Point2(0.5, 0.5);
		its.wi = wi;

		Float spp = m_size * m_size * m_sqrtSpp * m_sqrtSpp;
		Float normFactor = (double)m_size * m_size / spp * 0.25;
		
		for (int i = 0; i < m_size * m_sqrtSpp; i++) {
			for (int j = 0; j < m_size * m_sqrtSpp; j++) {
				BSDFSamplingRecord bRec(its, m_sampler);
				Spectrum spec = bsdf->sample(bRec, m_sampler->next2D());
				if (Frame::cosTheta(bRec.wo) <= 0)
					continue;

				int c = math::clamp(math::floorToInt((bRec.wo.x + 1.0) * 0.5 * m_size), 0, m_size - 1);
				int r = math::clamp(math::floorToInt((bRec.wo.y + 1.0) * 0.5 * m_size), 0, m_size - 1);

				bsdfValues[m_size - r - 1][c] += Vector3d(spec[0], spec[1], spec[2]) * 
					normFactor * Frame::cosTheta(bRec.wo);
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

MTS_EXPORT_UTILITY(ValidateSampleBSDF, "Visualize BSDF eval")
MTS_NAMESPACE_END