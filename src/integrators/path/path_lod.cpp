#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderproc.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

static StatsCounter avgPathLength("Path LoD tracer", "Average path length", EAverage);

class PathLodTracer : public MonteCarloIntegrator {
public:
	PathLodTracer(const Properties &props)
		: MonteCarloIntegrator(props) {
		std::string interpName = props.getString("interpType", "nearest");
		if (interpName == "nearest")
			m_interp = 0;
		else if (interpName == "linear")
			m_interp = 1;
		else if (interpName == "stoc_linear")
			m_interp = 2;
	}

	/// Unserialize from a binary data stream
	PathLodTracer(Stream *stream, InstanceManager *manager)
		: MonteCarloIntegrator(stream, manager) { }

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
		if (sampler->getClass()->getName() != "IndependentSampler")
			Log(EError, "The footprint sub-integrator should only be "
			"used in conjunction with the independent sampler");
		if (!m_subIntegrator->preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID))
			return false;

		Log(EInfo, "Start precompute pixel footprint");
		Scene *rScene = static_cast<Scene *>(Scheduler::getInstance()->getResource(sceneResID));
		Film *film = sensor->getFilm();
		film->clear();
		m_subIntegrator->render(rScene, queue, job, sceneResID,
			sensorResID, samplerResID);

		Log(EInfo, "Finish calling subIntegrator");

		m_footprint = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, film->getCropSize());
		film->develop(Point2i(0, 0), film->getSize(), Point2i(0, 0), m_footprint);

		float *data = m_footprint->getFloat32Data();
		Vector2i filmSize = film->getCropSize();
				
		ref<FileStream> stream = new FileStream("test_footprint.exr", FileStream::ETruncWrite);
		m_footprint->write(Bitmap::EOpenEXR, stream);

		return true;
	}

	bool render(Scene *scene,
		RenderQueue *queue, const RenderJob *job,
		int sceneResID, int sensorResID, int samplerResID) {
		ref<Scheduler> sched = Scheduler::getInstance();
		ref<Sensor> sensor = static_cast<Sensor *>(sched->getResource(sensorResID));
		ref<Film> film = sensor->getFilm();
		film->clear();

		m_lodLevel = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, film->getCropSize());

		size_t nCores = sched->getCoreCount();
		const Sampler *sampler = static_cast<const Sampler *>(sched->getResource(samplerResID, 0));
		size_t sampleCount = sampler->getSampleCount();

		int integratorResID = sched->registerResource(this);
		
		Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SIZE_T_FMT
			" %s, " SSE_STR ") ..", film->getCropSize().x, film->getCropSize().y,
			sampleCount, sampleCount == 1 ? "sample" : "samples", nCores,
			nCores == 1 ? "core" : "cores");

		/* This is a sampling-based integrator - parallelize */
		ref<ParallelProcess> proc = new BlockedRenderProcess(job,
			queue, scene->getBlockSize());
		
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

		ref<FileStream> stream = new FileStream("test_lod_level.exr", FileStream::ETruncWrite);
		m_lodLevel->write(Bitmap::EOpenEXR, stream);

		return proc->getReturnStatus() == ParallelProcess::ESuccess;
	}

	void renderBlockX(const Scene *scene,
		const Sensor *sensor, Sampler *sampler, ImageBlock *block,
		const bool &stop, const std::vector< TPoint2<uint8_t> > &points) {

		Float diffScaleFactor = 1.0f /
			std::sqrt((Float) sampler->getSampleCount());

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

		for (size_t i = 0; i<points.size(); ++i) {
			Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
			if (stop)
				break;

			sampler->generate(offset);

			// pixel footprint
			Spectrum fval = m_footprint->getPixel(offset);
			Float uvFootprint = fval[0];
			Float pixelChosenLevel = 0.0;

			if (m_interp == 0) {
				// nearest
				for (size_t j = 0; j<sampler->getSampleCount(); j++) {
					rRec.newQuery(queryType, sensor->getMedium());
					Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

					if (needsApertureSample)
						apertureSample = rRec.nextSample2D();
					if (needsTimeSample)
						timeSample = rRec.nextSample1D();

					Spectrum spec = sensor->sampleRayDifferential(
						sensorRay, samplePos, apertureSample, timeSample);

					sensorRay.scaleDifferential(diffScaleFactor);

					Float chosenLevel = 0;
					int levelOffset = 0;
					spec *= Li(sensorRay, rRec, uvFootprint, levelOffset, chosenLevel);
					block->put(samplePos, spec, rRec.alpha);
					sampler->advance();

					pixelChosenLevel += chosenLevel;
				}
			} else if (m_interp == 1) {
				// linear
				for (size_t j = 0; j<sampler->getSampleCount(); j++) {
					rRec.newQuery(queryType, sensor->getMedium());
					Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

					if (needsApertureSample)
						apertureSample = rRec.nextSample2D();
					if (needsTimeSample)
						timeSample = rRec.nextSample1D();

					Spectrum spec0 = sensor->sampleRayDifferential(
						sensorRay, samplePos, apertureSample, timeSample);

					sensorRay.scaleDifferential(diffScaleFactor);

					Float chosenLevel = 0;
					int levelOffset = 0;
					spec0 *= Li(sensorRay, rRec, uvFootprint, levelOffset, chosenLevel);

					Float w = chosenLevel - math::floorToInt(chosenLevel);
					spec0 *= (1.0 - w);

					// second lod
					rRec.newQuery(queryType, sensor->getMedium());
					Spectrum spec1 = sensor->sampleRayDifferential(
						sensorRay, samplePos, apertureSample, timeSample);

					sensorRay.scaleDifferential(diffScaleFactor);
					
					levelOffset = 1;
					spec1 *= Li(sensorRay, rRec, uvFootprint, levelOffset, chosenLevel);
					spec1 *= w;
					
					block->put(samplePos, spec0 + spec1, rRec.alpha);
					sampler->advance();

					pixelChosenLevel += chosenLevel;
				}
			} else if (m_interp == 2) {
				// stoc_linear
				Float w = 1.0;
				for (size_t j = 0; j<sampler->getSampleCount(); j++) {
					rRec.newQuery(queryType, sensor->getMedium());
					Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

					if (needsApertureSample)
						apertureSample = rRec.nextSample2D();
					if (needsTimeSample)
						timeSample = rRec.nextSample1D();

					Spectrum spec = sensor->sampleRayDifferential(
						sensorRay, samplePos, apertureSample, timeSample);

					sensorRay.scaleDifferential(diffScaleFactor);

					Float chosenLevel = 0;
					int levelOffset;

					if ((Float)j / sampler->getSampleCount() <= w)
						levelOffset = 0;
					else
						levelOffset = 1;

					spec *= Li(sensorRay, rRec, uvFootprint, levelOffset, chosenLevel);

					if (j < 10)
						w = (w * j + (1.0 + math::floorToInt(chosenLevel) - chosenLevel)) / (j + 1);

					block->put(samplePos, spec, rRec.alpha);
					sampler->advance();

					pixelChosenLevel += chosenLevel;
				}
			}

			m_lodLevel->setPixel(offset, Spectrum(pixelChosenLevel / (Float)sampler->getSampleCount()));
		}
	}

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec,
		Float uvFootprint, int levelOffset, Float &chosenLevel) {
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		RayDifferential ray(r);
		Spectrum Li(0.0f);
		bool scattered = false;

		/* Perform the first ray intersection (or ignore if the
		   intersection has already been provided). */
		ray.setUvFootprint(uvFootprint, levelOffset);
		rRec.rayIntersect(ray);
		ray.mint = Epsilon;

		Spectrum throughput(1.0f);
		Float eta = 1.0f;

		while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
			if (!its.isValid()) {
				/* If no intersection could be found, potentially return
				   radiance from a environment luminaire if it exists */
				if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
					&& (!m_hideEmitters || scattered))
					Li += throughput * scene->evalEnvironment(ray);
				break;
			}

			if (rRec.depth == 1) {
				chosenLevel = its.lodLevel;
			}

			const BSDF *bsdf = its.getBSDF(ray);

			/* Possibly include emitted radiance if requested */
			if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
				&& (!m_hideEmitters || scattered))
				Li += throughput * its.Le(-ray.d);

			/* Include radiance from a subsurface scattering model if requested */
			if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
				Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

			if ((rRec.depth >= m_maxDepth && m_maxDepth > 0)
				|| (m_strictNormals && dot(ray.d, its.geoFrame.n)
					* Frame::cosTheta(its.wi) >= 0)) {

				/* Only continue if:
				   1. The current path length is below the specifed maximum
				   2. If 'strictNormals'=true, when the geometric and shading
				      normals classify the incident direction to the same side */
				break;
			}

			/* ==================================================================== */
			/*                     Direct illumination sampling                     */
			/* ==================================================================== */

			/* Estimate the direct illumination if this is requested */
			DirectSamplingRecord dRec(its);
			dRec.setUvFootprint(uvFootprint, levelOffset);

			if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
				(bsdf->getType() & BSDF::ESmooth)) {
				Spectrum value = scene->sampleEmitterDirect(dRec, rRec.nextSample2D());
				if (!value.isZero()) {
					const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

					/* Allocate a record for querying the BSDF */
					BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

					/* Evaluate BSDF * cos(theta) */
					const Spectrum bsdfVal = bsdf->eval(bRec);

					/* Prevent light leaks due to the use of shading normals */
					if (!bsdfVal.isZero() && (!m_strictNormals
							|| dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {

						/* Calculate prob. of having generated that direction
						   using BSDF sampling */
						Float bsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
							? bsdf->pdf(bRec) : 0;

						/* Weight using the power heuristic */
						Float weight = miWeight(dRec.pdf, bsdfPdf);
						Li += throughput * value * bsdfVal * weight;
					}
				}
			}

			/* ==================================================================== */
			/*                            BSDF sampling                             */
			/* ==================================================================== */

			/* Sample BSDF * cos(theta) */
			Float bsdfPdf;
			BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
			Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
			if (bsdfWeight.isZero())
				break;

			scattered |= bRec.sampledType != BSDF::ENull;

			/* Prevent light leaks due to the use of shading normals */
			const Vector wo = its.toWorld(bRec.wo);
			Float woDotGeoN = dot(its.geoFrame.n, wo);
			if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
				break;

			bool hitEmitter = false;
			Spectrum value;

			/* Trace a ray in this direction */
			ray = Ray(its.p, wo, ray.time);
			ray.setUvFootprint(uvFootprint, levelOffset);

			if (scene->rayIntersect(ray, its)) {
				/* Intersected something - check if it was a luminaire */
				if (its.isEmitter()) {
					value = its.Le(-ray.d);
					dRec.setQuery(ray, its);
					hitEmitter = true;
				}
			} else {
				/* Intersected nothing -- perhaps there is an environment map? */
				const Emitter *env = scene->getEnvironmentEmitter();

				if (env) {
					if (m_hideEmitters && !scattered)
						break;

					value = env->evalEnvironment(ray);
					if (!env->fillDirectSamplingRecord(dRec, ray))
						break;
					hitEmitter = true;
				} else {
					break;
				}
			}

			/* Keep track of the throughput and relative
			   refractive index along the path */
			throughput *= bsdfWeight;
			eta *= bRec.eta;

			/* If a luminaire was hit, estimate the local illumination and
			   weight using the power heuristic */
			if (hitEmitter &&
				(rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
				/* Compute the prob. of generating that direction using the
				   implemented direct illumination sampling technique */
				const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
					scene->pdfEmitterDirect(dRec) : 0;
				Li += throughput * value * miWeight(bsdfPdf, lumPdf);
			}

			/* ==================================================================== */
			/*                         Indirect illumination                        */
			/* ==================================================================== */

			/* Set the recursive query type. Stop if no surface was hit by the
			   BSDF sample or if indirect illumination was not requested */
			if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
				break;
			rRec.type = RadianceQueryRecord::ERadianceNoEmission;

			if (rRec.depth++ >= m_rrDepth) {
				/* Russian roulette: try to keep path weights equal to one,
				   while accounting for the solid angle compression at refractive
				   index boundaries. Stop with at least some probability to avoid
				   getting stuck (e.g. due to total internal reflection) */

				Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
				if (rRec.nextSample1D() >= q)
					break;
				throughput /= q;
			}
		}

		/* Store statistics */
		avgPathLength.incrementBase();
		avgPathLength += rRec.depth;

		return Li;
	}

	void renderBlock(const Scene *scene,
		const Sensor *sensor, Sampler *sampler, ImageBlock *block,
		const bool &stop, const std::vector< TPoint2<uint8_t> > &points) const {
	}

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		return Spectrum(0.f);
	}

	inline Float miWeight(Float pdfA, Float pdfB) const {
		pdfA *= pdfA;
		pdfB *= pdfB;
		return pdfA / (pdfA + pdfB);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();

		if (cClass->derivesFrom(MTS_CLASS(Integrator))) {
			if (!cClass->derivesFrom(MTS_CLASS(SamplingIntegrator)))
				Log(EError, "The sub-integrator must be derived from the class SamplingIntegrator");
			m_subIntegrator = static_cast<SamplingIntegrator *>(child);
		} else {
			Integrator::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		MonteCarloIntegrator::serialize(stream, manager);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "PathLodTracer[" << endl
			<< "  maxDepth = " << m_maxDepth << "," << endl
			<< "  rrDepth = " << m_rrDepth << "," << endl
			<< "  strictNormals = " << m_strictNormals << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()

private:
	ref<SamplingIntegrator> m_subIntegrator;
	ref<Bitmap> m_footprint;
	ref<Bitmap> m_lodLevel;
	
	// 0: nearest, 1: linear, 2: stoc_linear?
	int m_interp;
};

MTS_IMPLEMENT_CLASS_S(PathLodTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(PathLodTracer, "Path tracer for LoD");
MTS_NAMESPACE_END
