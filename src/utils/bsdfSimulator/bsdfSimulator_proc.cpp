#include <mitsuba/core/fstream.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include "bsdfSimulator_proc.h"

MTS_NAMESPACE_BEGIN

// WorkProcessor
BSDFRayTracer::BSDFRayTracer(const Vector &wi, int sqrtNumParticles, int size, 
		const AABB2 &aabb, int minDepth, int maxDepth, int shadowOption, 
		int useFullSphere, Float distGiTexelScale)
	: m_wi(wi), m_sqrtNumParticles(sqrtNumParticles), m_size(size), m_minDepth(minDepth), m_maxDepth(maxDepth), 
	m_aabb(aabb), m_shadowOption(shadowOption), m_useFullSphere(useFullSphere), m_distGiTexelScale(distGiTexelScale) { }

ref<WorkUnit> BSDFRayTracer::createWorkUnit() const {
	return new RangeWorkUnit();
}

ref<WorkResult> BSDFRayTracer::createWorkResult() const {
	return new MultiLobeDistribution(m_minDepth + 2, m_size, m_useFullSphere);
}

ref<WorkProcessor> BSDFRayTracer::clone() const {
	return new BSDFRayTracer(m_wi, m_sqrtNumParticles, m_size, m_aabb, m_minDepth, m_maxDepth, 
		m_shadowOption, m_useFullSphere, m_distGiTexelScale);
}

void BSDFRayTracer::prepare() {
	Scene *scene = static_cast<Scene *>(getResource("scene"));
	m_scene = new Scene(scene);
	m_sampler = static_cast<Sampler *>(getResource("sampler"));
	m_scene->setSampler(m_sampler);
}

void BSDFRayTracer::process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop) {
	const RangeWorkUnit *range = static_cast<const RangeWorkUnit *>(workUnit);
	MultiLobeDistribution *res = static_cast<MultiLobeDistribution *>(workResult);

#if defined(USE_SQUARE_CONCENTRIC)
	double normFactor = (double)m_size * m_size / (2.0 * M_PI * (double)m_sqrtNumParticles * m_sqrtNumParticles);
#else
	double normFactor = (double)m_size * m_size / ((double)m_sqrtNumParticles * (double)m_sqrtNumParticles) * 0.25;
#endif
	//m_sampler->generate(Point2i(0));

	res->clear();

	m_distGiRange.x = (m_aabb.max.x - m_aabb.min.x) * m_distGiTexelScale;
	m_distGiRange.y = (m_aabb.max.y - m_aabb.min.y) * m_distGiTexelScale;
	
	//Log(EInfo, "process start %d: (%d, %d)", Thread::getID(), range->getRangeStart(), range->getRangeEnd());
	for (size_t i = range->getRangeStart(); i <= range->getRangeEnd() && !stop; i++) {
		//Log(EInfo, "working on particle %d", i);

		double weight;
		Point o = sampleRayOrigin(i, weight);
		if (weight < 1e-5) {
			//Log(EInfo, "sample origin fail...");
			continue;
		}

		RayDifferential ray(o, -m_wi, 0);
		RadianceQueryRecord rRec(m_scene, m_sampler);
		rRec.type = RadianceQueryRecord::ERadiance;

		Intersection its;
		Spectrum throughput = sampleReflectance(ray, rRec, its);

// 		if (throughput[0] < -Epsilon) {
// 			//Log(EInfo, "invalid throughput");
// 			continue;
// 		}

		//Assert(ray.d.z > -Epsilon);
		if (ray.d.z < 0 && !m_useFullSphere)
			throughput = Spectrum(0.f);

		if (!throughput.isZero() && ray.d.z > 0) {
			if (m_shadowOption == 1) {
				if (its.isValid() && m_aabb.contains(Point2(its.p.x, its.p.y)))
					throughput = Spectrum(0.f);
			}
			else if (m_shadowOption == 2) {
				if (its.isValid())
					throughput = Spectrum(0.f);
			}
		}

		if (rRec.depth > 0) {
			if (rRec.depth <= m_minDepth)
				res->getLobe(rRec.depth - 1)->put(ray.d, throughput, weight, normFactor);
			else
				res->getLobe(m_minDepth)->put(ray.d, throughput, weight, normFactor);
		}
		res->getLobe(m_minDepth + 1)->put(ray.d, throughput, weight, normFactor);

	}
	//Log(EInfo, "process done %d: (%d, %d)", Thread::getID(), range->getRangeStart(), range->getRangeEnd());
}

Point BSDFRayTracer::sampleRayOrigin(int idx, double &weight) {
	int r = idx / m_sqrtNumParticles;
	int c = idx % m_sqrtNumParticles;
	double x = m_aabb.min.x + (c + m_sampler->next1D()) / (double)m_sqrtNumParticles * (m_aabb.max.x - m_aabb.min.x);
	double y = m_aabb.min.y + (r + m_sampler->next1D()) / (double)m_sqrtNumParticles * (m_aabb.max.y - m_aabb.min.y);

	Point o(x, y, 1e2);
	Ray ray(o, Vector(0, 0, -1.0f), 0);
	Intersection its;
	m_scene->rayIntersect(ray, its);
	Normal normal = its.shFrame.n;

	double cosG = std::max(0.0, normal.z);
	if (cosG < 1e-5) {
		weight = 0.0;
		//Log(EInfo, "bad cosG: %.6f, %.6f, %.6f", normal.x, normal.y, normal.z);
		return o;
	}

	weight = std::max(0.0, dot(normal, m_wi)) / cosG;
	if (weight > 0.0) {
		ray = Ray(its.p + m_wi * ShadowEpsilon, m_wi, 0);
		o = its.p + m_wi * ShadowEpsilon * 10;

		// by default: global masking
		if (m_wi.z > 0) {
			if (m_scene->rayIntersect(ray)) {
				weight = 0.0;
				//Log(EInfo, "masked...");
			}
		} else {
			m_scene->rayIntersect(ray, its);
			
			if (its.isValid() && std::abs(its.p.x - x) < m_distGiRange.x && std::abs(its.p.y - y) < m_distGiRange.y) {
				weight = 0.0;
				//Log(EInfo, "%.6f, %.6f", std::abs(its.p.x - x) / (m_aabb.max.x - m_aabb.min.x) * 4096,
				//	std::abs(its.p.y - y) / (m_aabb.max.y - m_aabb.min.y) * 4096);
			}
		}

		/*
		if (m_shadowOption == 1) {
			m_scene->rayIntersect(ray, its);
			if (its.isValid() && m_aabb.contains(Point2(its.p.x, its.p.y))) {
				weight = 0.0;
				//Log(EInfo, "masked...");
			}
		}
		else if (m_shadowOption == 2) {
			if (m_scene->rayIntersect(ray)) {
				weight = 0.0;
				//Log(EInfo, "masked...");
			}
		}
		*/
	}

	return o;

// 	Point o(x, y, 10.f);
// 	o += m_wi * 1e3;
// 	Ray ray(o, -m_wi, 0);
// 	success = m_scene->rayIntersect(ray);
	
	//Intersection its;
	//success = m_scene->rayIntersect(ray, its);
	//Log(EInfo, "%d", m_aabb.contains(Point2(x, y)));
	//Log(EInfo, "(%.6f, %.6f), (%.6f, %.6f, %.6f)", x, y,
	//	its.p.x, its.p.y, its.p.z);

	//return o;
}

Spectrum BSDFRayTracer::sampleReflectance(RayDifferential &ray, RadianceQueryRecord &rRec, Intersection &getIts) {
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

			if (!its.isValid() || !m_aabb.contains(Point2(its.p.x, its.p.y)) ||
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

			/*
			if (rRec.depth > 0) {
				Vector wiMacro = its.baseFrame.toLocal(-ray.d);
				if (wiMacro.z > 0) {
					Log(EInfo, "Macro = %.6f, %.6f, %.6f", wiMacro.x, wiMacro.y, wiMacro.z);
					Log(EInfo, "Local = %.6f, %.6f, %.6f", its.wi.x, its.wi.y, its.wi.z);
				}
			}
			*/

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

BSDFRayTracer::BSDFRayTracer(Stream *stream, InstanceManager *manager) {
	for (int i = 0; i < 3; i++)
		m_wi[i] = stream->readFloat();
	m_sqrtNumParticles = stream->readInt();
	m_size = stream->readInt();
	m_minDepth = stream->readInt();
	m_maxDepth = stream->readInt();
	m_aabb = AABB2(stream);
	m_shadowOption = stream->readInt();
	m_useFullSphere = stream->readInt();
	m_distGiTexelScale = stream->readFloat();
}

void BSDFRayTracer::serialize(Stream *stream, InstanceManager *manager) const {
	for (int i = 0; i < 3; i++)
		stream->writeFloat(m_wi[i]);
	stream->writeInt(m_sqrtNumParticles);
	stream->writeInt(m_size);
	stream->writeInt(m_minDepth);
	stream->writeInt(m_maxDepth);
	m_aabb.serialize(stream);
	stream->writeInt(m_shadowOption);
	stream->writeInt(m_useFullSphere);
	stream->writeFloat(m_distGiTexelScale);
}

// ParallelProcess
BSDFSimulatorProcess::BSDFSimulatorProcess(const Vector &wi, int sqrtNumParticles,
		int size, const AABB2 &aabb, int minDepth, int maxDepth, int shadowOption, 
		int useFullSphere, Float distGiTexelScale)
		: m_wi(wi), m_sqrtNumParticles(sqrtNumParticles), m_size(size), m_minDepth(minDepth),
		m_maxDepth(maxDepth), m_aabb(aabb), m_shadowOption(shadowOption), m_useFullSphere(useFullSphere),
		m_distGiTexelScale(distGiTexelScale) {
	m_numParticles = m_sqrtNumParticles * m_sqrtNumParticles;
	m_start = 0;
	m_granularity = std::max((size_t)1, m_numParticles 
		/ (16 * Scheduler::getInstance()->getWorkerCount()));
	
	m_res = new MultiLobeDistribution(m_minDepth + 2, m_size, m_useFullSphere);
	m_res->clear();

	m_resultCount = 0;
	m_progress = new ProgressReporter("Rendering", (m_numParticles - 1) / m_granularity + 1, NULL);
	m_resultMutex = new Mutex();
}

BSDFSimulatorProcess::~BSDFSimulatorProcess() {
	if (m_progress)
		delete m_progress;
}

void BSDFSimulatorProcess::setGranularity(int granularity) {
	m_granularity = granularity;
}

ParallelProcess::EStatus BSDFSimulatorProcess::generateWork(WorkUnit *unit, int worker) {
	RangeWorkUnit &range = *static_cast<RangeWorkUnit *>(unit);

	if (m_start >= m_numParticles)
		return EFailure;

	int ed = std::min(m_start + m_granularity - 1, m_numParticles - 1);
	range.setRange(m_start, ed);
	m_start = ed + 1;

	return ESuccess;
}

void BSDFSimulatorProcess::processResult(const WorkResult *result, bool cancelled) {
	if (cancelled)
		return;

	const MultiLobeDistribution *res = static_cast<const MultiLobeDistribution *>(result);
	UniqueLock lock(m_resultMutex);

	m_res->put(res);
	m_progress->update(++m_resultCount);

	lock.unlock();
}

ref<WorkProcessor> BSDFSimulatorProcess::createWorkProcessor() const {
	return new BSDFRayTracer(m_wi, m_sqrtNumParticles, m_size, m_aabb, m_minDepth, m_maxDepth, m_shadowOption,
		m_useFullSphere, m_distGiTexelScale);
}

MTS_IMPLEMENT_CLASS_S(BSDFRayTracer, false, WorkProcessor)
MTS_IMPLEMENT_CLASS(BSDFSimulatorProcess, false, ParallelProcess)
MTS_NAMESPACE_END
