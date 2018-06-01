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
		m_wi = Vector(std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]));
		m_sqrtNumParticles = std::atoi(argv[5]);
		m_size = std::atoi(argv[6]);	
		
		m_xmin = std::atof(argv[7]);
		m_xmax = std::atof(argv[8]);
		m_ymin = std::atof(argv[9]);
		m_ymax = std::atof(argv[10]);

		// Example: minDepth = 2, it will create 4 lobes:
		// 1st-order, 2nd-order, (3 to maxDepth)-th order, and all orders
		m_minDepth = std::atoi(argv[11]);
		m_maxDepth = std::atoi(argv[12]);
		m_shadowOption = std::atoi(argv[13]);

		if (argc > 15)
			m_useFullSphere = std::atoi(argv[15]);
		else
			m_useFullSphere = 0;

		if (argc > 16)
			m_distGiTexelScale = std::atof(argv[16]);

		ParameterMap params;
		/*
		if (argc > 16) {
			params["ssR"] = argv[16];
			params["ssG"] = argv[17];
			params["ssB"] = argv[18];
		}
		if (argc > 19) {
			params["giR"] = argv[19];
			params["giG"] = argv[20];
			params["giB"] = argv[21];
		}
		if (argc > 22)
			params["ssExp"] = argv[22];
		if (argc > 23)
			params["giExp"] = argv[23];
		*/
		m_scene = loadScene(argv[1], params);
		
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
			m_minDepth, m_maxDepth, m_shadowOption, m_useFullSphere, m_distGiTexelScale);
		proc->bindResource("scene", sceneResID);
		proc->bindResource("sampler", samplerResID);
		m_scene->bindUsedResources(proc);

		Log(EInfo, "Start rendering.");

		sched->schedule(proc);
		sched->wait(proc);

		Log(EInfo, "Finish rendering.");

		// Normalization issues!
		double totValidParticles = (double)proc->m_res->getLobe(m_minDepth + 1)->m_totValidParticles;
		//double totValidParticles = m_numParticles;
		double totParticles = m_numParticles;
		double pScale = totParticles / totValidParticles;
		double totWeight = proc->m_res->getLobe(m_minDepth + 1)->m_totWeight;

		int numLobes = m_minDepth + 2;

		char txtFilename[256];
		FILE *fp;

// 		sprintf(txtFilename, "%s.txt", argv[14]);
// 		fp = fopen(txtFilename, "w");

		for (int i = 0; i < numLobes; i++) {
			char filename[256];
			if (m_minDepth == m_maxDepth && i + 1 != m_minDepth)
				continue;

			if (i < numLobes - 1)
				sprintf(filename, "%s_order_%d.exr", argv[14], i + 1);
			else
				sprintf(filename, "%s_order_all.exr", argv[14]);

			SphericalDistribution *lobe = proc->m_res->getLobe(i);
			
			//proc->m_res->getLobe(i)->scale(pScale);
			
			Float invTotWeight = 0.0;
			if (totWeight > 1e-8)
				invTotWeight = 1.0 / totWeight;
			lobe->scale(totParticles * invTotWeight);
			lobe->saveExr(fs::path(filename));

			//lobe->m_totValue /= totValidParticles;
			lobe->m_totValue *= invTotWeight;

			Vector3d &totalThroughput = lobe->m_totValue;
			Log(EInfo, "Total valid particles = %d / (%.0f, %d)", lobe->m_totValidParticles, 
				totValidParticles, m_numParticles);
			Log(EInfo, "Total thr = (%.6f, %.6f, %.6f)", totalThroughput[0], totalThroughput[1], totalThroughput[2]);
			
			/*fprintf(fp, "%.6f %.6f %.6f\n", totalThroughput[0], totalThroughput[1], totalThroughput[2]);*/
		}		
		/*fclose(fp);*/

// 		sprintf(txtFilename, "%s_moments.txt", argv[14]);
// 		fp = fopen(txtFilename, "w");
// 		for (int i = 0; i < 3; i++) {
// 			Vector3d moment = proc->m_res->getLobe(numLobes - 1)->m_moments[i];
// 			moment /= totValidParticles;
// 			Log(EInfo, "%d-th moment = (%.6f, %.6f, %.6f)", i, moment[0], moment[1], moment[2]);
// 
// 			fprintf(fp, "%.6f %.6f %.6f\n", moment[0], moment[1], moment[2]);
// 		}
// 		fclose(fp);
		
		//validate();
		return 0;
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

	int m_useFullSphere;
	Float m_distGiTexelScale;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(BSDFSimulator, "BSDF simulator")
MTS_NAMESPACE_END
