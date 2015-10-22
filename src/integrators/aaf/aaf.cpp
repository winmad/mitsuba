/*
	Added by Lifan Wu
	Oct 22, 2015
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>
#include <boost/filesystem/path.hpp>

MTS_NAMESPACE_BEGIN

class AAFIntegrator : public SamplingIntegrator {
public:
	AAFIntegrator(const Properties &props) : SamplingIntegrator(props) {
		m_initSamples = props.getSize("initSamples", 16);
		m_verbose = props.getBoolean("verbose", false);
	}

	AAFIntegrator(Stream *stream, InstanceManager *manager)
		: SamplingIntegrator(stream, manager) {
		m_subIntegrator = static_cast<SamplingIntegrator *>(manager->getInstance(stream));
		m_initSamples = stream->readInt();
		m_verbose = false;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();

		if (cClass->derivesFrom(MTS_CLASS(Integrator))) {
			if (!cClass->derivesFrom(MTS_CLASS(SamplingIntegrator)))
				Log(EError, "The sub-integrator must be derived from the class SamplingIntegrator");
			m_subIntegrator = static_cast<SamplingIntegrator *>(child);
		}
		else {
			Integrator::addChild(name, child);
		}
	}

	void configureSampler(const Scene *scene, Sampler *sampler) {
		SamplingIntegrator::configureSampler(scene, sampler);
		m_subIntegrator->configureSampler(scene, sampler);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int sensorResID, int samplerResID) {
		if (!SamplingIntegrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID))
			return false;
		if (m_subIntegrator == NULL)
			Log(EError, "No sub-integrator was specified!");
		Sampler *sampler = static_cast<Sampler *>(Scheduler::getInstance()->getResource(samplerResID, 0));
		Sensor *sensor = static_cast<Sensor *>(Scheduler::getInstance()->getResource(sensorResID));
// 		if (sampler->getClass()->getName() != "IndependentSampler")
// 			Log(EError, "The error-controlling integrator should only be "
// 			"used in conjunction with the independent sampler");
		if (!m_subIntegrator->preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID))
			return false;

		Vector2i filmSize = sensor->getFilm()->getSize();
		bool needsApertureSample = sensor->needsApertureSample();
		bool needsTimeSample = sensor->needsTimeSample();
		const int nSamples = 0;
		Float luminance = 0;

		Point2 apertureSample(0.5f);
		Float timeSample = 0.5f;
		RadianceQueryRecord rRec(scene, sampler);

		/* Estimate the overall luminance on the image plane */
		for (int i = 0; i<nSamples; ++i) {
			sampler->generate(Point2i(0));

			rRec.newQuery(RadianceQueryRecord::ERadiance, sensor->getMedium());
			rRec.extra = RadianceQueryRecord::EAdaptiveQuery;

			Point2 samplePos(rRec.nextSample2D());
			samplePos.x *= filmSize.x;
			samplePos.y *= filmSize.y;

			if (needsApertureSample)
				apertureSample = rRec.nextSample2D();
			if (needsTimeSample)
				timeSample = rRec.nextSample1D();

			RayDifferential eyeRay;
			Spectrum sampleValue = sensor->sampleRay(
				eyeRay, samplePos, apertureSample, timeSample);

			sampleValue *= m_subIntegrator->Li(eyeRay, rRec);
			luminance += sampleValue.getLuminance();
		}

		return true;
	}

	void renderBlock(const Scene *scene, const Sensor *sensor,
		Sampler *sampler, ImageBlock *block, const bool &stop,
		const std::vector< TPoint2<uint8_t> > &points) const {
		Float diffScaleFactor = 1.0f /
			std::sqrt((Float)sampler->getSampleCount());

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

		for (size_t i = 0; i < points.size(); ++i) {
			Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
			if (stop)
				break;

			sampler->generate(offset);

			for (size_t j = 0; j < sampler->getSampleCount(); j++) {
				rRec.newQuery(queryType, sensor->getMedium());
				Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

				if (needsApertureSample)
					apertureSample = rRec.nextSample2D();
				if (needsTimeSample)
					timeSample = rRec.nextSample1D();

				Spectrum spec = sensor->sampleRayDifferential(
					sensorRay, samplePos, apertureSample, timeSample);

				sensorRay.scaleDifferential(diffScaleFactor);

				spec *= m_subIntegrator->Li(sensorRay, rRec);
				block->put(samplePos, spec, rRec.alpha);
				sampler->advance();
			}
		}
	}

	Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
		return m_subIntegrator->Li(ray, rRec);
	}

	Spectrum E(const Scene *scene, const Intersection &its, const Medium *medium,
		Sampler *sampler, int nSamples, bool includeIndirect) const {
		return m_subIntegrator->E(scene, its, medium,
			sampler, nSamples, includeIndirect);
	}

	void postprocess(const Scene *scene, RenderQueue *queue,
		const RenderJob *job, int sceneResID, int sensorResID,
		int samplerResID) {
		const Film *film = scene->getSensor()->getFilm();
		const Bitmap *img = film->getImageBlock()->getBitmap();
		//boost::filesystem::path output{ "test.exr" };
		//Log(EInfo, "output path = %s", output.string().c_str());
		//img->write(output);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		SamplingIntegrator::serialize(stream, manager);
		manager->serialize(stream, m_subIntegrator.get());
		
		stream->writeInt(m_initSamples);
	}

	void bindUsedResources(ParallelProcess *proc) const {
		m_subIntegrator->bindUsedResources(proc);
	}

	void wakeup(ConfigurableObject *parent,
		std::map<std::string, SerializableObject *> &params) {
		m_subIntegrator->wakeup(this, params);
	}

	void cancel() {
		SamplingIntegrator::cancel();
		m_subIntegrator->cancel();
	}

	const Integrator *getSubIntegrator(int idx) const {
		if (idx != 0)
			return NULL;
		return m_subIntegrator.get();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "AdaptiveIntegrator[" << endl
			<< "  initSamples = " << m_initSamples << "," << endl
			<< "  subIntegrator = " << indent(m_subIntegrator->toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<SamplingIntegrator> m_subIntegrator;
	int m_initSamples;

	// 0: d1 (distance from light to pixel), 1: d2_min, 2: d2_max (distance from light to occluders)
	std::vector<Spectrum> info;

	bool m_verbose;
};

MTS_IMPLEMENT_CLASS_S(AAFIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(AAFIntegrator, "Axis-aligned filtering based integrator");
MTS_NAMESPACE_END
