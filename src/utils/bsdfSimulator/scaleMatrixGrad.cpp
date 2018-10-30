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
#include "scaleMatrixGrad_proc.h"

MTS_NAMESPACE_BEGIN

class ScaleMatrixGrad : public Utility {
public:
	int run(int argc, char **argv) {
		char *sceneFilepath = argv[1];
		m_wi = Vector(std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]));
		m_sqrtNumParticles = std::atoi(argv[5]);
		m_wiResolution = std::atoi(argv[6]);
		m_woResolution = std::atoi(argv[7]);
		
		m_xmin = std::atof(argv[8]);
		m_xmax = std::atof(argv[9]);
		m_ymin = std::atof(argv[10]);
		m_ymax = std::atof(argv[11]);

		m_maxDepth = std::atoi(argv[12]);
		m_shadowOption = std::atoi(argv[13]);

		char *weightFilepath = argv[14];
		char *outputPrefix = argv[15];

		if (argc > 16)
			m_distGiTexelScale = std::atof(argv[16]);

		m_scene = loadScene(fs::path(sceneFilepath));
		m_weightMap = new Bitmap(fs::path(weightFilepath));

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

		ref<ScaleMatrixGradProcess> proc = new ScaleMatrixGradProcess(m_wi,
			m_sqrtNumParticles, m_wiResolution, m_woResolution, AABB2(Point2(m_xmin, m_ymin), Point2(m_xmax, m_ymax)), 
			m_maxDepth, m_shadowOption, m_weightMap.get(), m_distGiTexelScale);
		proc->bindResource("scene", sceneResID);
		proc->bindResource("sampler", samplerResID);
		m_scene->bindUsedResources(proc);

		Log(EInfo, "Start rendering.");

		sched->schedule(proc);
		sched->wait(proc);

		Log(EInfo, "Finish rendering.");

		// Normalization issues!
		double totWeight = proc->m_res->getWeight();
		Float invTotWeight = 0.0;
		if (totWeight > 1e-8)
			invTotWeight = 1.0 / totWeight;

		char filename[256];
		sprintf(filename, "%s.exr", outputPrefix);
		const std::vector<Float> &grad = proc->m_res->getData();
		outputGrad(grad, m_numParticles * invTotWeight, filename);
		
		int totElems = 2 * m_wiResolution * m_wiResolution * 2 * m_woResolution * m_woResolution;
		Vector3 totalThroughput(grad[3 * totElems], grad[3 * totElems + 1], grad[3 * totElems + 2]);
		totalThroughput *= (double)m_numParticles * invTotWeight;
		Log(EInfo, "Total thr = (%.6f, %.6f, %.6f)", totalThroughput[0], totalThroughput[1], totalThroughput[2]);

		return 0;		
	}

	void outputGrad(const std::vector<Float> &grad, Float scale, char *filename) {
		ref<Bitmap> bitmap = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, 
			Vector2i(2 * m_wiResolution * m_wiResolution, 2 * m_woResolution * m_woResolution + 1));
		float *data = bitmap->getFloat32Data();

		int matHeight = 2 * m_wiResolution * m_wiResolution;
		int matWidth = 2 * m_woResolution * m_woResolution;
		for (int r = 0; r < matHeight; r++) {
			for (int c = 0; c < matWidth; c++) {
				for (int k = 0; k < 3; k++) {
					*data++ = grad[((r * matWidth) + c) * 3 + k] * scale;
				}
			}
		}

		for (int k = 0; k < 3; k++) {
			*data++ = grad[matWidth * matHeight * 3 + k] * scale;
		}

		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		bitmap->write(Bitmap::EOpenEXR, stream);
	}

	ref<Scene> m_scene;
	ref<Sampler> m_sampler;
	Vector m_wi;
	int m_numParticles, m_sqrtNumParticles;
	int m_wiResolution, m_woResolution;
	int m_minDepth;
	int m_maxDepth;
	double m_xmin, m_xmax, m_ymin, m_ymax;
	int m_shadowOption;

	ref<Bitmap> m_weightMap;
	Float m_distGiTexelScale;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(ScaleMatrixGrad, "BSDF simulator")
MTS_NAMESPACE_END
