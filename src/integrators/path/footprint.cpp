#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderproc.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

class FootprintIntegrator : public SamplingIntegrator {
public:
	FootprintIntegrator(const Properties &props)
		: SamplingIntegrator(props) { }

	/// Unserialize from a binary data stream
	FootprintIntegrator(Stream *stream, InstanceManager *manager)
		: SamplingIntegrator(stream, manager) { }

	bool render(Scene *scene,
		RenderQueue *queue, const RenderJob *job,
		int sceneResID, int sensorResID, int samplerResID) {
		ref<Scheduler> sched = Scheduler::getInstance();
		ref<Sensor> sensor = static_cast<Sensor *>(sched->getResource(sensorResID));
		ref<Film> film = sensor->getFilm();
		film->clear();

		size_t nCores = sched->getCoreCount();
		const Sampler *sampler = static_cast<const Sampler *>(sched->getResource(samplerResID, 0));
		size_t sampleCount = 4;

		Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SIZE_T_FMT
			" %s, " SSE_STR ") ..", film->getCropSize().x, film->getCropSize().y,
			sampleCount, sampleCount == 1 ? "sample" : "samples", nCores,
			nCores == 1 ? "core" : "cores");
		
		/* This is a sampling-based integrator - parallelize */
		ref<ParallelProcess> proc = new BlockedRenderProcess(job,
			queue, scene->getBlockSize());
		int integratorResID = sched->registerResource(this);		
		proc->bindResource("integrator", integratorResID);
		proc->bindResource("scene", sceneResID);
		proc->bindResource("sensor", sensorResID);
		proc->bindResource("sampler", samplerResID);
		scene->bindUsedResources(proc);
		bindUsedResources(proc);
		sched->schedule(proc);

		m_process = proc;
		sched->wait(proc);
		m_process = NULL;
		sched->unregisterResource(integratorResID);

		return proc->getReturnStatus() == ParallelProcess::ESuccess;
	}

	void renderBlock(const Scene *scene,
		const Sensor *sensor, Sampler *sampler, ImageBlock *block,
		const bool &stop, const std::vector< TPoint2<uint8_t> > &points) const {

		//Float diffScaleFactor = 1.0f /
		//	std::sqrt((Float) sampler->getSampleCount());
		Float diffScaleFactor = 1.0f;

		bool needsApertureSample = sensor->needsApertureSample();
		bool needsTimeSample = sensor->needsTimeSample();

		RadianceQueryRecord rRec(scene, sampler);
		Point2 apertureSample(0.5f);
		Float timeSample = 0.5f;
		RayDifferential sensorRay;

		block->clear();

		uint32_t queryType = RadianceQueryRecord::ESensorRay;

		if (!sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we don't have to */
			queryType &= ~RadianceQueryRecord::EOpacity;

		int sampleCount = 4;
		std::vector<const Shape*> itsShapes(4);

		for (size_t i = 0; i<points.size(); ++i) {
			Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
			if (stop)
				break;

			sampler->generate(offset);
			
			Spectrum pixelRes(0.0);
			
			for (size_t j = 0; j<sampleCount; j++) {
				rRec.newQuery(queryType, sensor->getMedium());
				Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

				if (needsApertureSample)
					apertureSample = rRec.nextSample2D();
				if (needsTimeSample)
					timeSample = rRec.nextSample1D();

				Spectrum spec = sensor->sampleRayDifferential(
					sensorRay, samplePos, apertureSample, timeSample);

				sensorRay.scaleDifferential(diffScaleFactor);
				sensorRay.setUvFootprint(1.f, 0);

				Spectrum li = Li(sensorRay, rRec);
				itsShapes[j] = rRec.its.shape;
				
				pixelRes += li;
				//block->put(samplePos, spec, rRec.alpha);
				sampler->advance();

				//Log(EInfo, "(%d, %d), %.6f, %.6f, %.6f", offset.x, offset.y,
				//	spec[0], spec[1], spec[2]);
			}

			bool atBoundary = false;
			for (int j = 0; j < sampleCount && !atBoundary; j++) {
				if (itsShapes[j] == NULL) {
					atBoundary = true;
					break;
				}
				for (int k = j + 1; k < sampleCount; k++) {
					if (itsShapes[j] != itsShapes[k]) {
						atBoundary = true;
						break;
					}
				}
			}

			if (atBoundary) {
				pixelRes = Spectrum(0.0);
			}
			block->put(Point2(offset) + Vector2(0.5), pixelRes / sampleCount, rRec.alpha);
		}
	}

	Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
		Spectrum result(0.0);

		if (!rRec.rayIntersect(ray))
			return result;

		Intersection &its = rRec.its;

		its.computePartials(ray);
		result[0] = std::min(sqrt(its.dudx * its.dudx + its.dvdx * its.dvdx),
			sqrt(its.dudy * its.dudy + its.dvdy * its.dvdy));
		
		result[1] = result[0] * its.effTexReso;
		result[2] = its.lodLevel;

		return result;
	}

	void renderBlockX(const Scene *scene,
		const Sensor *sensor, Sampler *sampler, ImageBlock *block,
		const bool &stop, const std::vector< TPoint2<uint8_t> > &points) {
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		SamplingIntegrator::serialize(stream, manager);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "FootprintIntegrator[" << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()

private:
	ref<Bitmap> m_footprint;
};

MTS_IMPLEMENT_CLASS_S(FootprintIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(FootprintIntegrator, "Path tracer for LoD");
MTS_NAMESPACE_END
