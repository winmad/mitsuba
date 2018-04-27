#include <mitsuba/core/fstream.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include "lobeDerivative_proc.h"

MTS_NAMESPACE_BEGIN

// WorkProcessor
BSDFDerivativeRayTracer::BSDFDerivativeRayTracer(int numVars, const Vector &wi, int sqrtNumParticles, int size,
		const AABB2 &aabb, int minDepth, int maxDepth, int shadowOption)
	: m_numVars(numVars), m_wi(wi), m_sqrtNumParticles(sqrtNumParticles), 
	m_size(size), m_minDepth(minDepth), m_maxDepth(maxDepth),
	m_aabb(aabb), m_shadowOption(shadowOption) { }

void BSDFDerivativeRayTracer::sampleDerivative(RayDifferential &ray,
		RadianceQueryRecord &rRec, Intersection &getIts, 
		Float normFactor, MultiLobeDistribution *res) {
	const Scene *scene = rRec.scene;
	Intersection &its = rRec.its;

	rRec.rayIntersect(ray);
	Spectrum throughput(1.0f);

	Spectrum sumDiff[4];
	for (int i = 0; i < 4; i++) {
		sumDiff[i] = Spectrum(0.0f);
	}

	while (rRec.depth <= m_maxDepth) {
		getIts = its;
		if (throughput.isZero())
			break;

		if (!its.isValid() || !m_aabb.contains(Point2(its.p.x, its.p.y)) ||
			rRec.depth == m_maxDepth) {
			if (ray.d.z < 0)
				break;
			if (m_shadowOption == 1) {
				if (its.isValid() && m_aabb.contains(Point2(its.p.x, its.p.y)))
					break;
			}
			else if (m_shadowOption == 2) {
				if (its.isValid())
					break;
			}

			// put result
			for (int i = 0; i < 4; i++) {
				Spectrum tmp = throughput * sumDiff[i];
				res->getLobe(i)->put(ray.d, tmp, 1.0, normFactor);
			}
			return;
		}

		const BSDF *bsdf = its.getBSDF();
		BSDFSamplingRecord bRec(its, rRec.sampler);

		Spectrum bsdfVal = bsdf->sample(bRec, rRec.nextSample2D());
		throughput *= bsdfVal;
		if (bsdfVal.isZero())
			break;

		// update sumDiff
		for (int i = 0; i < 2; i++) {
			if (!bRec.used[i])
				continue;
			for (int c = 0; c < 3; c++) {
				sumDiff[i][c] += bRec.dAlbedo[i][c];
				sumDiff[i + 2][c] += bRec.dRoughness[i];
			}
		}

		const Vector wo = its.toWorld(bRec.wo);
		ray = Ray(its.p, wo, ray.time);

		scene->rayIntersect(ray, its);
		
		rRec.depth++;
	}
}

void BSDFDerivativeRayTracer::process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop) {
	const RangeWorkUnit *range = static_cast<const RangeWorkUnit *>(workUnit);
	MultiLobeDistribution *res = static_cast<MultiLobeDistribution *>(workResult);

	double normFactor = (double)m_size * m_size / ((double)m_sqrtNumParticles * (double)m_sqrtNumParticles) * 0.25;
	//m_sampler->generate(Point2i(0));

	res->clear();

	//Log(EInfo, "process start %d: (%d, %d)", Thread::getID(), range->getRangeStart(), range->getRangeEnd());
	for (size_t i = range->getRangeStart(); i <= range->getRangeEnd() && !stop; i++) {
		//Log(EInfo, "working on particle %d", i);

		bool success;
		Point o = sampleRayOrigin(i, success);
		if (!success) {
			//Log(EInfo, "sample origin fail...");
			continue;
		}

		RayDifferential ray(o, -m_wi, 0);
		RadianceQueryRecord rRec(m_scene, m_sampler);
		rRec.type = RadianceQueryRecord::ERadiance;

		Intersection its;
		sampleDerivative(ray, rRec, its, 
			normFactor, res);
	}
	//Log(EInfo, "process done %d: (%d, %d)", Thread::getID(), range->getRangeStart(), range->getRangeEnd());
}

Point BSDFDerivativeRayTracer::sampleRayOrigin(int idx, bool &success) {
	int r = idx / m_sqrtNumParticles;
	int c = idx % m_sqrtNumParticles;
	double x = m_aabb.min.x + (c + m_sampler->next1D()) / (double)m_sqrtNumParticles * (m_aabb.max.x - m_aabb.min.x);
	double y = m_aabb.min.y + (r + m_sampler->next1D()) / (double)m_sqrtNumParticles * (m_aabb.max.y - m_aabb.min.y);

	Point o(x, y, 10.0f);
	o += m_wi * 1e3;
	Ray ray(o, -m_wi, 0);
	success = m_scene->rayIntersect(ray);

	return o;
}

ref<WorkUnit> BSDFDerivativeRayTracer::createWorkUnit() const {
	return new RangeWorkUnit();
}

ref<WorkResult> BSDFDerivativeRayTracer::createWorkResult() const {
	return new MultiLobeDistribution(m_numVars, m_size);
}

ref<WorkProcessor> BSDFDerivativeRayTracer::clone() const {
	return new BSDFDerivativeRayTracer(m_numVars, m_wi, m_sqrtNumParticles,
		m_size, m_aabb, m_minDepth, m_maxDepth, m_shadowOption);
}

void BSDFDerivativeRayTracer::prepare() {
	Scene *scene = static_cast<Scene *>(getResource("scene"));
	m_scene = new Scene(scene);
	m_sampler = static_cast<Sampler *>(getResource("sampler"));
	m_scene->setSampler(m_sampler);
}

BSDFDerivativeRayTracer::BSDFDerivativeRayTracer(Stream *stream, InstanceManager *manager) {
	m_numVars = stream->readInt();
	for (int i = 0; i < 3; i++)
		m_wi[i] = stream->readFloat();
	m_sqrtNumParticles = stream->readInt();
	m_size = stream->readInt();
	m_aabb = AABB2(stream);
	m_minDepth = stream->readInt();
	m_maxDepth = stream->readInt();
	m_shadowOption = stream->readInt();
}

void BSDFDerivativeRayTracer::serialize(Stream *stream, InstanceManager *manager) const {
	stream->writeInt(m_numVars);
	for (int i = 0; i < 3; i++)
		stream->writeFloat(m_wi[i]);
	stream->writeInt(m_sqrtNumParticles);
	stream->writeInt(m_size);
	m_aabb.serialize(stream);
	stream->writeInt(m_minDepth);
	stream->writeInt(m_maxDepth);
	stream->writeInt(m_shadowOption);
}

// ParallelProcess
BSDFDerivativeProcess::BSDFDerivativeProcess(int numVars, const Vector &wi, int sqrtNumParticles,
		int size, const AABB2 &aabb, int minDepth, int maxDepth, int shadowOption)
		: m_numVars(numVars), m_wi(wi), m_sqrtNumParticles(sqrtNumParticles), m_size(size), 
		m_minDepth(minDepth), m_maxDepth(maxDepth), m_aabb(aabb), m_shadowOption(shadowOption) {
	m_numParticles = m_sqrtNumParticles * m_sqrtNumParticles;
	m_start = 0;
	m_granularity = std::max((size_t)1, m_numParticles
		/ (16 * Scheduler::getInstance()->getWorkerCount()));

	m_res = new MultiLobeDistribution(m_numVars, m_size);
	m_res->clear();

	m_resultCount = 0;
	m_progress = new ProgressReporter("Rendering", (m_numParticles - 1) / m_granularity + 1, NULL);
	m_resultMutex = new Mutex();
}

BSDFDerivativeProcess::~BSDFDerivativeProcess() {
	if (m_progress)
		delete m_progress;
}

void BSDFDerivativeProcess::setGranularity(int granularity) {
	m_granularity = granularity;
}

ParallelProcess::EStatus BSDFDerivativeProcess::generateWork(WorkUnit *unit, int worker) {
	RangeWorkUnit &range = *static_cast<RangeWorkUnit *>(unit);

	if (m_start >= m_numParticles)
		return EFailure;

	int ed = std::min(m_start + m_granularity - 1, m_numParticles - 1);
	range.setRange(m_start, ed);
	m_start = ed + 1;

	return ESuccess;
}

void BSDFDerivativeProcess::processResult(const WorkResult *result, bool cancelled) {
	if (cancelled)
		return;

	const MultiLobeDistribution *res = static_cast<const MultiLobeDistribution *>(result);

	m_res->put(res);

	m_progress->update(++m_resultCount);
	UniqueLock lock(m_resultMutex);
	lock.unlock();
}

ref<WorkProcessor> BSDFDerivativeProcess::createWorkProcessor() const {
	return new BSDFDerivativeRayTracer(m_numVars, m_wi, m_sqrtNumParticles, m_size,
		m_aabb, m_minDepth, m_maxDepth, m_shadowOption);
}

MTS_IMPLEMENT_CLASS_S(BSDFDerivativeRayTracer, false, WorkProcessor)
MTS_IMPLEMENT_CLASS(BSDFDerivativeProcess, false, ParallelProcess)
MTS_NAMESPACE_END
