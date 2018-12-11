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

class EffBrdfEstimator : public Utility {
public:
	int run(int argc, char **argv) {
		char *inputParamFilename = argv[2];
		m_sqrtSpp = std::atof(argv[3]);

		// Example: minDepth = 2, it will create 2 lobes:
		// (1 to minDepth)-th order, (1 to maxDepth)-th order
		m_minDepth = std::atoi(argv[4]);
		m_maxDepth = std::atoi(argv[5]);
		m_shadowOption = std::atoi(argv[6]);
		m_maskingOption = std::atoi(argv[7]);

		// for GI neighborhood
		m_distGiTexelScale = std::atof(argv[8]);
		m_outputPrefix = std::string(argv[9]);

		std::vector<AABB2> aabbs;
		std::vector<Vector> wis;
		std::vector<Vector> wos;

		FILE* fp = fopen(inputParamFilename, "r");
		double inputVal[10];
		while (fscanf(fp, "%lf ", &inputVal[0]) != EOF) {
			for (int i = 1; i < 10; i++)
				fscanf(fp, "%lf ", &inputVal[i]);
			aabbs.push_back(AABB2(Point2(inputVal[0], inputVal[2]), Point2(inputVal[1], inputVal[3])));
			wis.push_back(Vector(inputVal[4], inputVal[5], inputVal[6]));
			wos.push_back(Vector(inputVal[7], inputVal[8], inputVal[9]));
		}
		fclose(fp);

		/*
		m_xmin = std::atof(argv[2]);
		m_xmax = std::atof(argv[3]);
		m_ymin = std::atof(argv[4]);
		m_ymax = std::atof(argv[5]);

		m_wi = Vector(std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]));
		m_wo = Vector(std::atof(argv[9]), std::atof(argv[10]), std::atof(argv[11]));

		m_sqrtSpp = std::atof(argv[12]);

		// Example: minDepth = 2, it will create 2 lobes:
		// (1 to minDepth)-th order, (1 to maxDepth)-th order
		m_minDepth = std::atoi(argv[13]);
		m_maxDepth = std::atoi(argv[14]);
		m_shadowOption = std::atoi(argv[15]);
		m_maskingOption = std::atoi(argv[16]);

		// for GI neighborhood
		m_distGiTexelScale = std::atof(argv[17]);

		m_outputPrefix = std::string(argv[18]);
		*/

		ParameterMap params;
		if (argc > 10) {
			params["angular"] = argv[10];
			params["spatial"] = argv[11];
		}
		m_scene = loadScene(argv[1], params);
		m_scene->initialize();

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

		char filename[256];
		Log(EInfo, "Start working");
		std::vector<Spectrum> results;

		for (int i = 0; i < aabbs.size(); i++) {
			// init
			m_aabb = aabbs[i];
			m_wi = wis[i];
			m_wo = wos[i];
			
			intEffBrdfOverP(m_visProjArea, m_effBrdf);

			// output eff brdf
			Float totArea = 0.f;
			Spectrum res(0.f);
			for (int i = 0; i < m_sqrtSpp; i++) {
				for (int j = 0; j < m_sqrtSpp; j++) {
					res += m_effBrdf[i][j];
					totArea += m_visProjArea[i][j];
				}
			}
			if (totArea > 1e-4) {
				res /= totArea;
			}

			results.push_back(res);

// 			Log(EInfo, "Total_vis_area = %.6f", totArea / (m_sqrtSpp * m_sqrtSpp));
// 			Log(EInfo, "Eff_BRDF_value = (%.6f, %.6f, %.6f)", res[0], res[1], res[2]);

// 			sprintf(filename, "%s_value.txt", m_outputPrefix.c_str());
// 			fp = fopen(filename, "w");
// 			fprintf(fp, "%.8f %.8f %.8f", res[0], res[1], res[2]);
// 			fclose(fp);
// 
// 			sprintf(filename, "%s_order_all.exr", m_outputPrefix.c_str());
// 			outputEffBrdf(m_effBrdf, filename);
		}

		Log(EInfo, "Done!");

		sprintf(filename, "%s_values.txt", m_outputPrefix.c_str());
		fp = fopen(filename, "w");
		for (int i = 0; i < results.size(); i++) {
			fprintf(fp, "%.8f %.8f %.8f\n", results[i][0], results[i][1], results[i][2]);
		}
		fclose(fp);
		
		return 0;
	}

	void intEffBrdfOverP(std::vector<std::vector<Float> > &visProjArea,
		std::vector<std::vector<Spectrum> > &res) {
		Vector2 distGiRange;
		distGiRange.x = (m_aabb.max.x - m_aabb.min.x) * m_distGiTexelScale;
		distGiRange.y = (m_aabb.max.y - m_aabb.min.y) * m_distGiTexelScale;

		res.resize(m_sqrtSpp);
		for (int i = 0; i < m_sqrtSpp; i++)
			res[i].resize(m_sqrtSpp, Spectrum(0.f));

		visProjArea.resize(m_sqrtSpp);
		for (int i = 0; i < m_sqrtSpp; i++)
			visProjArea[i].resize(m_sqrtSpp, 0.f);

#pragma omp parallel for
		for (int s = 0; s < m_sqrtSpp * m_sqrtSpp; s++) {
			ref<Sampler> sampler = m_samplers[Thread::getID() % 233];

			int i = s / m_sqrtSpp;
			int j = s % m_sqrtSpp;

			Float w;
			RayDifferential ray = samplePathStart(i, j, m_aabb, sampler, w);

			if (w < 1e-5) {
				// masked
				continue;
			}

			RadianceQueryRecord rRec(m_scene, sampler);
			rRec.type = RadianceQueryRecord::ERadiance;
			Intersection its;

			AABB2 giRangeAABB(Point2(ray.o.x - distGiRange.x, ray.o.y - distGiRange.y),
				Point2(ray.o.x + distGiRange.x, ray.o.y + distGiRange.y));
			Spectrum effBrdf = sampleReflectance(giRangeAABB, ray, rRec, its);

			res[i][j] = effBrdf * w;
			visProjArea[i][j] = w;
		}
	}

	Spectrum sampleReflectance(const AABB2 &giRangeAABB, 
		RayDifferential &ray, RadianceQueryRecord &rRec, Intersection &getIts) {
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		MediumSamplingRecord mRec;

		rRec.rayIntersect(ray);
		Spectrum throughput(1.0f);
		Spectrum res(0.0f);

		while (rRec.depth <= m_maxDepth) {
			getIts = its;
			if (throughput.isZero())
				break;

			if (!its.isValid() || !giRangeAABB.contains(Point2(its.p.x, its.p.y)) ||
				rRec.depth == m_maxDepth) {
// 				if (rRec.depth == 0)
// 				if (true)
// 					Log(EInfo, "%d, (%.6f, %.6f, %.6f), %d", rRec.depth, its.p.x, its.p.y, its.p.z,
// 						giRangeAABB.contains(Point2(its.p.x, its.p.y)));
				return res;
			}

			const BSDF *bsdf = its.getBSDF();

			// Next event estimation
			bool shadowed = false;
			Ray shadowRay(its.p + m_wo * ShadowEpsilon, m_wo, 0);
			if (m_shadowOption == 1) {
				Intersection shadowIts;
				scene->rayIntersect(shadowRay, shadowIts);
				if (shadowIts.isValid() && m_aabb.contains(Point2(shadowIts.p.x, shadowIts.p.y)))
					shadowed = true;
			} else if (m_shadowOption == 2) {
				if (scene->rayIntersect(shadowRay))
					shadowed = true;
			}

			if (!shadowed) {
				BSDFSamplingRecord bRecEval(its, its.toLocal(m_wo));
				Spectrum bsdfEval = bsdf->eval(bRecEval);
				res += throughput * bsdfEval;
// 				Log(EInfo, "thr = (%.6f, %.6f, %.6f)", throughput[0], throughput[1], throughput[2]);
// 				Log(EInfo, "bsdfEval = (%.6f, %.6f, %.6f)", bsdfEval[0], bsdfEval[1], bsdfEval[2]);
			}
			
			// Sample next direction
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
			ray = Ray(its.p, wo, ray.time);

			scene->rayIntersect(ray, its);
			rRec.depth++;
		}

		return res;
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
		w = std::max(0.0, dot(normal, m_wi)) / cosG;
		if (w > 0.0) {
			ray = Ray(its.p + m_wi * ShadowEpsilon, m_wi, 0);
			o = its.p + m_wi * ShadowEpsilon * 100.0;

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

		return RayDifferential(o, -m_wi, 0);
	}

	void outputEffBrdf(std::vector<std::vector<Spectrum> > &res, char *filename) const {
		ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, Vector2i(m_sqrtSpp, m_sqrtSpp));
		float *data = bitmap->getFloat32Data();
		for (int r = 0; r < m_sqrtSpp; r++) {
			for (int c = 0; c < m_sqrtSpp; c++) {
				*data++ = res[r][c][0];
				*data++ = res[r][c][1];
				*data++ = res[r][c][2];
			}
		}
		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		bitmap->write(Bitmap::EOpenEXR, stream);
	}

	ref<Scene> m_scene;
	std::string m_outputPrefix;
	ref_vector<Sampler> m_samplers;
	int m_sqrtSpp;
	int m_minDepth;
	int m_maxDepth;
	double m_xmin, m_xmax, m_ymin, m_ymax;
	Vector3 m_wi, m_wo;
	AABB2 m_aabb;

	// 1: local, 2: global
	int m_shadowOption, m_maskingOption;

	Float m_distGiTexelScale;

	std::vector<std::vector<Spectrum> > m_effBrdf;
	std::vector<std::vector<Float> > m_visProjArea;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(EffBrdfEstimator, "Effective BRDF estimator")
MTS_NAMESPACE_END
