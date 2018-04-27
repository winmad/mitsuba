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
#include "lobeDerivative_proc.h"

MTS_NAMESPACE_BEGIN

class LobeDerivative : public Utility {
public:
	int run(int argc, char **argv) {
		// 0, 1: albedo
		// 2, 3: roughness
		m_numVars = 4;
		bool varMask[4] = { false, true, false, true };

		m_wi = Vector(std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]));
		m_sqrtNumParticles = std::atoi(argv[5]);
		m_size = std::atoi(argv[6]);

		m_xmin = std::atof(argv[7]);
		m_xmax = std::atof(argv[8]);
		m_ymin = std::atof(argv[9]);
		m_ymax = std::atof(argv[10]);

		m_minDepth = std::atoi(argv[11]);
		m_maxDepth = std::atoi(argv[12]);
		m_shadowOption = std::atoi(argv[13]);

		ParameterMap params;
		if (argc > 15) {
			params["ssR"] = argv[15];
			params["ssG"] = argv[16];
			params["ssB"] = argv[17];
		}
		if (argc > 18) {
			params["giR"] = argv[18];
			params["giG"] = argv[19];
			params["giB"] = argv[20];
		}
		if (argc > 21)
			params["ssExp"] = argv[21];
		if (argc > 22)
			params["giExp"] = argv[22];
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

		ref<BSDFDerivativeProcess> proc = new BSDFDerivativeProcess(
			m_numVars, m_wi, m_sqrtNumParticles, m_size, 
			AABB2(Point2(m_xmin, m_ymin), Point2(m_xmax, m_ymax)),
			m_minDepth, m_maxDepth, m_shadowOption);
		proc->bindResource("scene", sceneResID);
		proc->bindResource("sampler", samplerResID);
		m_scene->bindUsedResources(proc);

		Log(EInfo, "Start rendering.");

		sched->schedule(proc);
		sched->wait(proc);

		Log(EInfo, "Finish rendering.");

		//double totValidParticles = 0.0;
		//for (int i = 0; i <= m_minDepth; i++)
		//	totValidParticles += (double)proc->m_res->getLobe(i)->m_totValidParticles;
		double totValidParticles = m_numParticles;
		int numLobes = m_numVars;

		char txtFilename[256];
		FILE *fp;
		
		sprintf(txtFilename, "%s.txt", argv[14]);
		fp = fopen(txtFilename, "w");
		for (int i = 0; i < numLobes; i++) {
			if (!varMask[i])
				continue;
			char filename[256];
			sprintf(filename, "%s_var_%d.exr", argv[14], i);

			proc->m_res->getLobe(i)->saveExr(fs::path(filename));

			proc->m_res->getLobe(i)->m_totValue /= totValidParticles;
			Vector3d &totalThroughput = proc->m_res->getLobe(i)->m_totValue;
			Log(EInfo, "Total valid particles = %d / (%.0f, %d)", proc->m_res->getLobe(i)->m_totValidParticles,
				totValidParticles, m_numParticles);
			Log(EInfo, "Total derivative = (%.6f, %.6f, %.6f)", totalThroughput[0], totalThroughput[1], totalThroughput[2]);

			fprintf(fp, "%.6f %.6f %.6f\n", totalThroughput[0], totalThroughput[1], totalThroughput[2]);
		}
		fclose(fp);

// 		sprintf(txtFilename, "%s_moments_derivative.txt", argv[14]);
// 		fp = fopen(txtFilename, "w");
// 		for (int i = 0; i < 3; i++) {
// 			// i-th moment
// 			for (int k = 0; k < numLobes; k++) {
// 				if (!varMask[k])
// 					continue;
// 				// variable k
// 				Vector3d moment = proc->m_res->getLobe(k)->m_moments[i];
// 				moment /= totValidParticles;
// 				Log(EInfo, "%d-th moment, var %d = (%.6f, %.6f, %.6f)", i, k, moment[0], moment[1], moment[2]);
// 
// 				fprintf(fp, "%.6f %.6f %.6f ", moment[0], moment[1], moment[2]);
// 			}
// 			fprintf(fp, "\n");
// 		}
// 		fclose(fp);
	
		return 0;
	}
	
public:
	ref<Scene> m_scene;
	ref<Sampler> m_sampler;
	Vector m_wi;
	int m_numParticles, m_sqrtNumParticles;
	int m_size;
	int m_minDepth;
	int m_maxDepth;
	double m_xmin, m_xmax, m_ymin, m_ymax;
	int m_shadowOption;
	int m_numVars;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(LobeDerivative, "Lobe derivative estimator")
MTS_NAMESPACE_END
