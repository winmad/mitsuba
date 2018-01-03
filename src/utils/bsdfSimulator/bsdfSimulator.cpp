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
#include "bsdfSimulator_proc.h"

MTS_NAMESPACE_BEGIN

class BSDFSimulator : public Utility {
public:
	int run(int argc, char **argv) {
		m_scene = loadScene(argv[1]);		
		m_wi = Vector(std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]));
		m_sqrtNumParticles = std::atoi(argv[5]);
		m_size = std::atoi(argv[6]);	
		
		m_xmin = std::atof(argv[7]);
		m_xmax = std::atof(argv[8]);
		m_ymin = std::atof(argv[9]);
		m_ymax = std::atof(argv[10]);

		// Example: minDepth = 2, it will create 3 lobes:
		// 1st-order, 2nd-order, and (3 to maxDepth)-th order
		m_minDepth = std::atoi(argv[11]);
		m_maxDepth = std::atoi(argv[12]);
		m_shadowOption = std::atoi(argv[13]);
		
		// init
		m_scene->initialize();
		m_wi = normalize(m_wi);
		m_numParticles = m_sqrtNumParticles * m_sqrtNumParticles;

		Properties props = Properties("independent");
		props.setInteger("seed", 19931004);
		m_sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), props));
		m_sampler->configure();

		//runSingleThread(argv);

		ref<Scheduler> sched = Scheduler::getInstance();
		int sceneResID = sched->registerResource(m_scene);

		std::vector<SerializableObject *> samplers(sched->getCoreCount());
		for (size_t i = 0; i < sched->getCoreCount(); i++) {
			ref<Sampler> clonedSampler = m_sampler->clone();
			clonedSampler->incRef();
			samplers[i] = clonedSampler.get();
		}
		int samplerResID = sched->registerMultiResource(samplers);
		for (size_t i = 0; i < samplers.size(); i++)
			samplers[i]->decRef();

		ref<BSDFSimulatorProcess> proc = new BSDFSimulatorProcess(m_wi,
			m_sqrtNumParticles, m_size, AABB2(Point2(m_xmin, m_ymin), Point2(m_xmax, m_ymax)), 
			m_minDepth, m_maxDepth, m_shadowOption);
		proc->bindResource("scene", sceneResID);
		proc->bindResource("sampler", samplerResID);
		m_scene->bindUsedResources(proc);

		sched->schedule(proc);
		sched->wait(proc);

		Log(EInfo, "Finish rendering.");

		char txtFilename[256];
		memcpy(txtFilename, argv[14], sizeof(char) * strlen(argv[14]));
		strcat(txtFilename, ".txt");
		FILE *fp = fopen(txtFilename, "w");

		double totValidParticles = 0.0;
		for (int i = 0; i <= m_minDepth; i++)
			totValidParticles += (double)proc->m_res->getLobe(i)->m_totValidParticles;

		for (int i = 0; i <= m_minDepth; i++) {
			char filename[256];
			sprintf(filename, "%s_order_%d.exr", argv[14], i + 1);
			proc->m_res->getLobe(i)->saveExr(fs::path(filename));

			proc->m_res->getLobe(i)->m_totValue /= totValidParticles;
			Vector3d &totalThroughput = proc->m_res->getLobe(i)->m_totValue;
			Log(EInfo, "Total valid particles = %d / (%.0f, %d)", proc->m_res->getLobe(i)->m_totValidParticles, 
				totValidParticles, m_numParticles);
			Log(EInfo, "Total thr = (%.6f, %.6f, %.6f)", totalThroughput[0], totalThroughput[1], totalThroughput[2]);
			
			fprintf(fp, "%.6f %.6f %.6f\n", totalThroughput[0], totalThroughput[1], totalThroughput[2]);
		}		
		fclose(fp);
		
		//validate();
		return 0;
	}

	// not sync with the multi-thread version
	void runSingleThread(char **argv) {
		std::vector<std::vector<Vector3d> > bsdfValues(m_size, std::vector<Vector3d>(m_size));
		for (int i = 0; i < m_size; i++)
			for (int j = 0; j < m_size; j++)
				bsdfValues[i][j] = Vector3d(0.0f);
		Vector3d totalThroughput(0.0f);
		int totalValidParticles = 0;

		// (1 / numParticles) / (area / numBins) = (1 / numParticles) / (4.0 / (m_size * m_size))
		double normFactor = (double)m_size * m_size / (double)m_numParticles * 0.25;

		for (int i = 0; i < m_numParticles; i++) {
			if ((i + 1) % 100000 == 0) {
				std::cout << "Running " << (i + 1) << "/" << m_numParticles << "\r" << std::flush;
			}

			Point o = sampleRayOrigin(i);
			RayDifferential ray(o, -m_wi, 0);
			RadianceQueryRecord rRec(m_scene, m_sampler);
			rRec.type = RadianceQueryRecord::ERadiance;

			Spectrum throughput = sample(ray, rRec);

// 			if (throughput.isZero()) {
// 				Log(EInfo, "thr = (%.8f, %.8f, %.8f),\ndir = (%.6f, %.6f, %.6f)",
// 					throughput[0], throughput[1], throughput[2],
// 					ray.d.x, ray.d.y, ray.d.z);
// 			}

			if (throughput[0] < -Epsilon)
				continue;
			//Assert(ray.d.z > -Epsilon);
			if (ray.d.z < 0)
				throughput = Spectrum(0.f);
				//continue;

			totalValidParticles++;		

			int c = (ray.d.x > 0.9999f ? m_size - 1 : floor((ray.d.x + 1.0) * 0.5 * m_size));
			int r = (ray.d.y > 0.9999f ? m_size - 1 : floor((ray.d.y + 1.0) * 0.5 * m_size));

// 			Spectrum rgb;
// 			throughput.toLinearRGB(rgb[0], rgb[1], rgb[2]);
// 			Log(EInfo, "(%.6f, %.6f, %.6f), (%.6f, %.6f, %.6f)",
// 				throughput[0], throughput[1], throughput[2],
// 				rgb[0], rgb[1], rgb[2]);

			for (int k = 0; k < 3; k++) {
				bsdfValues[m_size - r - 1][c][k] += throughput[k] * ray.d.z * normFactor;
				//bsdfValues[m_size - r - 1][c][k] += throughput[k] * normFactor;
				totalThroughput[k] += throughput[k];
			}
		}
		totalThroughput /= (double)totalValidParticles;

		Log(EInfo, "Total valid particles = %d / %d", totalValidParticles, m_numParticles);
		Log(EInfo, "Total thr = (%.6f, %.6f, %.6f)", totalThroughput[0], totalThroughput[1], totalThroughput[2]);

		normFactor = (double)m_numParticles / (double)totalValidParticles;

		ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, Vector2i(m_size));
		float *data = bitmap->getFloat32Data();
		for (int r = 0; r < m_size; r++) {
			for (int c = 0; c < m_size; c++) {
				*data++ = bsdfValues[r][c][0] * normFactor;
				*data++ = bsdfValues[r][c][1] * normFactor;
				*data++ = bsdfValues[r][c][2] * normFactor;
			}
		}
		fs::path filename(argv[1]);
		filename.replace_extension(".exr");
		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		bitmap->write(Bitmap::EOpenEXR, stream);
	}

	Point sampleRayOrigin(int idx) {
		int r = idx / m_sqrtNumParticles;
		int c = idx % m_sqrtNumParticles;
		double x = m_xmin + (c + m_sampler->next1D()) / (double)m_sqrtNumParticles * (m_xmax - m_xmin);
		double y = m_ymin + (r + m_sampler->next1D()) / (double)m_sqrtNumParticles * (m_ymax - m_ymin);
		//double x = m_xmin + (c + 0.5) / (double)m_sqrtNumParticles * (m_xmax - m_xmin);
		//double y = m_ymin + (r + 0.5) / (double)m_sqrtNumParticles * (m_ymax - m_ymin);
		Point o(x, y, 0);
		o = o + m_wi * 1e3;
		return o;
	}

	Spectrum sample(RayDifferential &ray, RadianceQueryRecord &rRec) {
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		MediumSamplingRecord mRec;

		rRec.rayIntersect(ray);
		Spectrum throughput(1.0f);

		while (rRec.depth <= m_maxDepth) {
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

				if (!its.isValid()) {
					if (rRec.depth < m_minDepth)
						throughput = Spectrum(0.0f);
					return throughput;
				}

				//Log(EInfo, "=== depth %d ===", rRec.depth);
				//Log(EInfo, "(%.6f, %.6f, %.6f)", its.p.x, its.p.y, its.p.z);

				const BSDF *bsdf = its.getBSDF();
				BSDFSamplingRecord bRec(its, rRec.sampler);

				// Mark back-faced intersection as invalid
				if (Frame::cosTheta(bRec.wi) <= 0) {
					return Spectrum(-1.0f);
				}

				Spectrum bsdfVal = bsdf->sample(bRec, rRec.nextSample2D());
				if (bsdfVal.isZero()) {
// 					Log(EInfo, "zero bsdf term");
// 					Log(EInfo, "%d", rRec.depth);
// 					Log(EInfo, "cos(wi) = %.6f", Frame::cosTheta(bRec.wi));
// 					Log(EInfo, "fresnel = %.6f", fresnelConductorExact(Frame::cosTheta(bRec.wi),
// 						Spectrum(0.0f / 1.000277), Spectrum(1.0f / 1.000277)));
// 					Log(EInfo, "pos = (%.6f, %.6f, %.6f), n = (%.6f, %.6f, %.6f)",
// 						its.p.x, its.p.y, its.p.z,
// 						its.shFrame.n.x, its.shFrame.n.y, its.shFrame.n.z);
					break;
				}
				throughput *= bsdfVal;

				const Vector wo = its.toWorld(bRec.wo);
				if (its.isMediumTransition())
					rRec.medium = its.getTargetMedium(wo);
				ray = Ray(its.p, wo, ray.time);
				scene->rayIntersect(ray, its);
			}
			rRec.depth++;
		}
		return Spectrum(0.f);
	}

	void validate() {
		std::vector<std::vector<Vector3d> > bsdfValues(m_size, std::vector<Vector3d>(m_size));
		for (int i = 0; i < m_size; i++)
			for (int j = 0; j < m_size; j++)
				bsdfValues[i][j] = Vector3d(0.0f);

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

		Properties props = Properties("roughconductor");
		props.setString("distribution", "beckmann");
		props.setFloat("alpha", 1.0);
		props.setString("material", "none");
		BSDF *bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props));
		bsdf->configure();

		int sampleCount = 100;
#pragma omp parallel for
		for (int i = 0; i < m_size; i++) {
			for (int j = 0; j < m_size; j++) {
				Vector3d res(0.0);
				for (int k = 0; k < sampleCount; k++) {
					double x = (j + m_sampler->next1D()) / (double)m_size * 2.0 - 1.0;
					double y = (i + m_sampler->next1D()) / (double)m_size * 2.0 - 1.0;
					if (x * x + y * y > 1.0)
						continue;
					double z = sqrt(1.0 - x * x - y * y);

					Intersection its;
					BSDFSamplingRecord bRec(its, m_wi, Vector(x, y, z));
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
		fs::path filename("valid.exr");
		filename.replace_extension(".exr");
		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		bitmap->write(Bitmap::EOpenEXR, stream);
	}

	ref<Scene> m_scene;
	ref<Sampler> m_sampler;
	Vector m_wi;
	int m_numParticles, m_sqrtNumParticles;
	int m_size;
	int m_minDepth;
	int m_maxDepth;
	double m_xmin, m_xmax, m_ymin, m_ymax;
	int m_shadowOption;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(BSDFSimulator, "BSDF simulator")
MTS_NAMESPACE_END
