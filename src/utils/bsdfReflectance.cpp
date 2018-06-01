#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/aabb.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <boost/filesystem/path.hpp>

MTS_NAMESPACE_BEGIN

class BSDFReflectance : public Utility {
public:
	int run(int argc, char **argv) {
		m_xmin = std::atof(argv[2]);
		m_xmax = std::atof(argv[3]);
		m_ymin = std::atof(argv[4]);
		m_ymax = std::atof(argv[5]);
		m_sqrtSpp = std::atof(argv[6]);

		// Example: minDepth = 2, it will create 2 lobes:
		// (1 to minDepth)-th order, (1 to maxDepth)-th order
		m_minDepth = std::atoi(argv[7]);
		m_maxDepth = std::atoi(argv[8]);
		m_shadowOption = std::atoi(argv[9]);
		m_maskingOption = std::atoi(argv[10]);

		// for GI neighborhood
		m_extendRange.x = std::atof(argv[11]);
		m_extendRange.y = std::atof(argv[12]);

		// spatially-varying blocks
		m_size.x = std::atoi(argv[13]);
		m_size.y = std::atoi(argv[14]);
		m_blockNums.x = std::atoi(argv[15]);
		m_blockNums.y = std::atoi(argv[16]);
		m_stIdx = std::atoi(argv[17]);
		m_edIdx = std::atoi(argv[18]);

 		ParameterMap params;
// 		if (argc > 15) {
// 			params["ssR"] = argv[15];
// 			params["ssG"] = argv[16];
// 			params["ssB"] = argv[17];
// 		}
// 		if (argc > 18) {
// 			params["giR"] = argv[18];
// 			params["giG"] = argv[19];
// 			params["giB"] = argv[20];
// 		}
// 		if (argc > 21)
// 			params["ssExp"] = argv[21];
// 		if (argc > 22)
// 			params["giExp"] = argv[22];
		m_scene = loadScene(argv[1], params);
		
		// init
		m_aabb = AABB2(Point2(m_xmin, m_ymin), Point2(m_xmax, m_ymax));
		m_scene->initialize();
		m_spp = m_sqrtSpp * m_sqrtSpp;

		// init samplers
		Properties props = Properties("independent");
		props.setInteger("seed", 19931004);
		m_samplers.resize(233);
		m_samplers[0] = static_cast<Sampler *>(PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));
		m_samplers[0]->configure();
		for (int i = 1; i < 233; i++) {
			m_samplers[i] = m_samplers[0]->clone();
		}

		m_rhoSingle.resize(m_blockNums.y);
		m_rhoAll.resize(m_blockNums.y);
		for (int i = 0; i < m_blockNums.y; i++) {
			m_rhoSingle[i].resize(m_blockNums.x, Spectrum(0.0));
			m_rhoAll[i].resize(m_blockNums.x, Spectrum(0.0));
		}

		Float xstep = (m_xmax - m_xmin) / m_blockNums.x;
		Float ystep = (m_ymax - m_ymin) / m_blockNums.y;

#pragma omp parallel for
		for (int r = 0; r < m_blockNums.y; r++) {
			ref<Sampler> sampler = m_samplers[Thread::getID() % 233];
			for (int c = 0; c < m_blockNums.x; c++) {
				int blockIdx = r * m_blockNums.x + c;
				if (blockIdx < m_stIdx || blockIdx >= m_edIdx)
					continue;

				Log(EInfo, "working on block (%d, %d)", c, r);

				AABB2 aabb(Point2(m_xmin + c * xstep, m_ymin + r * ystep),
					Point2(m_xmin + (c + 1) * xstep, m_ymax + (r + 1) * ystep));
				AABB2 aabbLimit(aabb);
				aabbLimit.min -= m_extendRange;
				aabbLimit.max += m_extendRange;

				// estimate meso-normal
				Vector nMeso(0.0);
				Vector2i spcnt;
				spcnt.x = math::clamp(m_size.x / m_blockNums.x, 16, 128);
				spcnt.y = math::clamp(m_size.y / m_blockNums.y, 16, 128);

				for (int i = 0; i < spcnt.y; i++) {
					for (int j = 0; j < spcnt.x; j++) {
						Float x = aabb.min.x + (j + sampler->next1D()) / spcnt.x * (aabb.max.x - aabb.min.x);
						Float y = aabb.min.y + (i + sampler->next1D()) / spcnt.y * (aabb.max.y - aabb.min.y);
						Point o(x, y, 1e2);
						Ray ray(o, Vector(0, 0, -1.0), 0);
						Intersection its;
						m_scene->rayIntersect(ray, its);
						nMeso += its.shFrame.n / its.shFrame.n.z;
					}
				}
				Float len = nMeso.length();
				if (len > 1e-5f)
					nMeso /= len;

				// Monte Carlo integration
				Spectrum &rhoSingle = m_rhoSingle[r][c];
				Spectrum &rhoAll = m_rhoAll[r][c];
				for (int i = 0; i < m_sqrtSpp; i++) {
					for (int j = 0; j < m_sqrtSpp; j++) {
						// uniformly sample position, sample wi by cos(wi)/pi
						Float w;
						RayDifferential ray = samplePathStart(i, j, aabb, sampler, w);

						if (w < 1e-5) {
							// masked
							continue;
						}

						Float cosWiMeso = dot(-ray.d, nMeso);
						if (cosWiMeso < 1e-5) {
							continue;
						}
						w *= nMeso.z / cosWiMeso * M_PI;

						RadianceQueryRecord rRec(m_scene, sampler);
						rRec.type = RadianceQueryRecord::ERadiance;
						Intersection its;
						Spectrum throughput = sampleReflectance(aabbLimit, ray, rRec, its);

						if (ray.d.z < 0)
							continue;

						if (!throughput.isZero()) {
							if (m_shadowOption == 1) {
								if (its.isValid() && aabbLimit.contains(Point2(its.p.x, its.p.y)))
									throughput = Spectrum(0.f);
							} else if (m_shadowOption == 2) {
								if (its.isValid())
									throughput = Spectrum(0.f);
							}
						}

						if (rRec.depth > 0 && rRec.depth <= m_minDepth) {
							rhoSingle += throughput * w;
						}
						rhoAll += throughput * w;
					}
				}
				rhoSingle /= m_spp;
				rhoAll /= m_spp;
			}
		}

		// output reflectance map
		char filename[256];
		sprintf(filename, "%s_rho_order_%d.exr", argv[19], m_minDepth);
		outputReflectanceMap(m_rhoSingle, filename);
		
		sprintf(filename, "%s_rho_order_all.exr", argv[19]);
		outputReflectanceMap(m_rhoAll, filename);

		//validate();
		return 0;
	}

	Spectrum sampleReflectance(AABB2 &aabbLimit, RayDifferential &ray, RadianceQueryRecord &rRec, Intersection &getIts) {
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		MediumSamplingRecord mRec;

		rRec.rayIntersect(ray);
		Spectrum throughput(1.0f);

		while (rRec.depth <= m_maxDepth) {
			getIts = its;
			if (throughput.isZero())
				break;

			if (rRec.medium && rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler)) {
				const PhaseFunction *phase = rRec.medium->getPhaseFunction();
				throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

				PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
				Float phaseVal = phase->sample(pRec, rRec.sampler);
				throughput *= phaseVal;
				if (phaseVal == 0)
					break;

				// Trace a ray
				ray = Ray(mRec.p, pRec.wo, ray.time);
				ray.mint = 0;
				scene->rayIntersect(ray, its);
			}
			else {
				if (rRec.medium)
					throughput *= mRec.transmittance / mRec.pdfFailure;

				if (!its.isValid() || !aabbLimit.contains(Point2(its.p.x, its.p.y)) ||
					rRec.depth == m_maxDepth) {
// 				if (rRec.depth == 0)
// 					Log(EInfo, "%d, (%.6f, %.6f, %.6f), %d", rRec.depth, its.p.x, its.p.y, its.p.z,
// 						m_aabb.contains(Point2(its.p.x, its.p.y)));
					return throughput;
				}

				const BSDF *bsdf = its.getBSDF();
				BSDFSamplingRecord bRec(its, rRec.sampler);

				// Mark back-faced intersection as invalid
				//if (Frame::cosTheta(bRec.wi) <= 0) {
				//	return Spectrum(-1.0f);
				//}

				Spectrum bsdfVal = bsdf->sample(bRec, rRec.nextSample2D());
				throughput *= bsdfVal;
				if (bsdfVal.isZero()) {
					//Log(EInfo, "zero bsdf, wo = (%.6f, %.6f, %.6f)", bRec.wo.x, bRec.wo.y, bRec.wo.z);
					break;
				}

				const Vector wo = its.toWorld(bRec.wo);
				if (its.isMediumTransition())
					rRec.medium = its.getTargetMedium(wo);
				ray = Ray(its.p, wo, ray.time);

				scene->rayIntersect(ray, its);
			}
			rRec.depth++;
		}

		return throughput;
	}

	RayDifferential samplePathStart(int i, int j, const AABB2 &aabb, Sampler *sampler, double &w) {
		// sample position
		Float x = aabb.min.x + (j + sampler->next1D()) / m_sqrtSpp * (aabb.max.x - aabb.min.x);
		Float y = aabb.min.y + (i + sampler->next1D()) / m_sqrtSpp * (aabb.max.y - aabb.min.y);

		Point o(x, y, 1e2);
		Ray ray(o, Vector(0, 0, -1.0f), 0);
		Intersection its;
		m_scene->rayIntersect(ray, its);
		Normal normal = its.shFrame.n;

		double cosG = std::max(0.0, normal.z);
		if (cosG < 1e-5) {
			w = 0;
			return RayDifferential();
		}

		// sample wi
		Vector wi = warp::squareToCosineHemisphere(sampler->next2D());
		w = std::max(0.0, dot(normal, wi)) / cosG;
		if (w > 0.0) {
			ray = Ray(its.p + wi * ShadowEpsilon, wi, 0);
			o = its.p + wi * ShadowEpsilon * 10.0;

			if (m_maskingOption == 1) {
				m_scene->rayIntersect(ray, its);
				if (its.isValid() && aabb.contains(Point2(its.p.x, its.p.y))) {
					w = 0.0;
				}
			} else if (m_maskingOption == 2) {
				if (m_scene->rayIntersect(ray)) {
					w = 0.0;
				}
			}
		}

		return RayDifferential(o, -wi, 0);
	}

	void validate() {
		/*
		Properties props = Properties("diffuse");
		props.setSpectrum("reflectance", Spectrum(0.8f));
		BSDF *bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props));
		bsdf->configure();
		*/

		/*
		Properties props = Properties("roughconductor");
		props.setString("distribution", "beckmann");
		props.setFloat("alpha", 0.1);
		props.setString("material", "Cu");
		BSDF *bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props));
		bsdf->configure();
		*/

		/*
		Properties props = Properties("roughconductor");
		props.setString("distribution", "beckmann");
		props.setFloat("alpha", 1.0);
		props.setString("material", "none");
		BSDF *bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props));
		bsdf->configure();
		*/

		BSDF *bsdf = m_scene->getShapes()[0]->getBSDF();

		int sampleCount = 10;
		Vector3d res(0.0);
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				ref<Sampler> sampler = m_samplers[Thread::getID() % 233];
				for (int k = 0; k < sampleCount; k++) {
					double x = (j + sampler->next1D()) / (double)m_sqrtSpp;
					double y = (i + sampler->next1D()) / (double)m_sqrtSpp;
					Vector wi = warp::squareToUniformHemisphereConcentric(Point2(x, y));

					Float ox = m_aabb.min.x + (j + sampler->next1D()) / m_sqrtSpp * (m_aabb.max.x - m_aabb.min.x);
					Float oy = m_aabb.min.y + (i + sampler->next1D()) / m_sqrtSpp * (m_aabb.max.y - m_aabb.min.y);
					Point o(ox, oy, 1e2);
					Ray ray(o, Vector(0, 0, -1.0), 0);
					Intersection its;
					m_scene->rayIntersect(ray, its);

					BSDFSamplingRecord bRec(its, wi);
					Spectrum tmp = bsdf->sample(bRec, sampler->next2D());

					for (int c = 0; c < 3; c++) {
						res[c] += tmp[c] * wi.z * 2.0 * M_PI;
					}
				}
			}
		}
		res /= m_spp * sampleCount;
		Log(EInfo, "rho = (%.6f, %.6f, %.6f)", res[0], res[1], res[2]);
	}

	void outputReflectanceMap(std::vector<std::vector<Spectrum> > &rho, char *filename) const {
		ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, m_blockNums);
		float *data = bitmap->getFloat32Data();
		for (int r = 0; r < m_blockNums.y; r++) {
			for (int c = 0; c < m_blockNums.x; c++) {
				*data++ = rho[r][c][0];
				*data++ = rho[r][c][1];
				*data++ = rho[r][c][2];
			}
		}
		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		bitmap->write(Bitmap::EOpenEXR, stream);
	}

	ref<Scene> m_scene;
	ref_vector<Sampler> m_samplers;
	int m_sqrtSpp, m_spp;
	int m_minDepth;
	int m_maxDepth;
	double m_xmin, m_xmax, m_ymin, m_ymax;
	AABB2 m_aabb;
	
	// 1: local, 2: global
	int m_shadowOption, m_maskingOption;
	Vector2 m_extendRange;

	Vector2i m_size;
	Vector2i m_blockNums;
	int m_stIdx, m_edIdx;

	std::vector<std::vector<Spectrum> > m_rhoSingle, m_rhoAll;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(BSDFReflectance, "BSDF reflectance estimators")
MTS_NAMESPACE_END
