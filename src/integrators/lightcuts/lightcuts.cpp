/*
	Added by Lifan Wu
	Oct 18, 2015
*/

#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

class LightcutsIntegrator : public SamplingIntegrator {
public:
	LightcutsIntegrator(const Properties &props) : SamplingIntegrator(props) {
		m_shadingSamples = props.getSize("shadingSamples", 1);
		m_rayLength = props.getFloat("rayLength", -1);
	}

	/// Unserialize from a binary data stream
	LightcutsIntegrator(Stream *stream, InstanceManager *manager)
		: SamplingIntegrator(stream, manager) {
			m_shadingSamples = stream->readSize();
			m_rayLength = stream->readFloat();
			configure();
	}
	
	/// Serialize to a binary data stream
	void serialize(Stream *stream, InstanceManager *manager) const {
		SamplingIntegrator::serialize(stream, manager);
		stream->writeSize(m_shadingSamples);
		stream->writeFloat(m_rayLength);
	}

	void configureSampler(const Scene *scene, Sampler *sampler) {
		SamplingIntegrator::configureSampler(scene, sampler);

		if (m_shadingSamples > 1)
			sampler->request2DArray(m_shadingSamples);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue,
		const RenderJob *job, int sceneResID, int sensorResID,
		int samplerResID) {
		Integrator::preprocess(scene, queue, job, sceneResID,
				sensorResID, samplerResID);
		if (m_rayLength < 0) {
			m_rayLength = scene->getAABB().getBSphere().radius * 0.5f;
			Log(EInfo, "Setting occlusion ray length to %f", m_rayLength);
		}
		return true;
	}

	Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
		/* Some aliases and local variables */
		Spectrum Li(0.0f);
		Point2 sample;

		/* Perform the first ray intersection (or ignore if the
		   intersection has already been provided). */
		if (!rRec.rayIntersect(ray)) {
			/* If no intersection could be found, possibly return
			   radiance from a background emitter */
			return Spectrum(1.0f);
		}

		/* Figure out how many shading samples to take, and where the
		   required random numbers should come from */
		Point2 *sampleArray = &sample;
		size_t numShadingSamples = m_shadingSamples;

		bool adaptiveQuery = (rRec.extra & RadianceQueryRecord::EAdaptiveQuery);

		if (numShadingSamples > 1 && rRec.depth == 1 && !adaptiveQuery) {
			sampleArray = rRec.sampler->next2DArray(numShadingSamples);
		} else {
			/* This integrator is used recursively by another integrator.
			   Be less accurate as this sample will not directly be observed. */
			numShadingSamples = 1;
			sample = rRec.nextSample2D();
		}

		const Intersection &its = rRec.its;
		for (size_t i=0; i<numShadingSamples; ++i) {
			Vector d = its.toWorld(warp::squareToCosineHemisphere(sampleArray[i]));

			Ray shadowRay(its.p, d, Epsilon, m_rayLength, ray.time);
			if (!rRec.scene->rayIntersect(shadowRay))
				Li += Spectrum(1.0f);
		}

		Li /= static_cast<Float>(numShadingSamples);

		return Li;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "AmbientOcclusionIntegrator[" << endl
			<< "  shadingSamples = " << m_shadingSamples << "," << endl
			<< "  rayLength = " << m_rayLength << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	size_t m_shadingSamples;
	Float m_rayLength;
};

MTS_IMPLEMENT_CLASS_S(LightcutsIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(LightcutsIntegrator, "Lightcuts integrator");
MTS_NAMESPACE_END