#include <mitsuba/core/fstream.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include "bsdfSimulator_proc.h"

MTS_NAMESPACE_BEGIN

// WorkResult
SphericalDistribution::SphericalDistribution(int size) : m_size(size) {
	m_values = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat64, Vector2i(m_size));
}

void SphericalDistribution::clear() {
	m_totValue = Vector3d(0.0);
	m_totValidParticles = 0;
	double *data = m_values->getFloat64Data();
	for (int i = 0; i < m_size * m_size * SPECTRUM_SAMPLES; i++)
		*data++ = 0.0;
}

void SphericalDistribution::put(const SphericalDistribution *dist) {
	Assert(m_size == dist->m_size);
	m_totValue += dist->m_totValue;
	m_totValidParticles += dist->m_totValidParticles;
	double *target = m_values->getFloat64Data();
	const double *source = dist->m_values->getFloat64Data();
	for (int i = 0; i < m_size * m_size * SPECTRUM_SAMPLES; i++)
		*target++ += *source++;
}

void SphericalDistribution::put(const Vector &dir, const Spectrum &value, double normFactor) {
	int c = (dir.x > 0.9999f ? m_size - 1 : floor((dir.x + 1.0) * 0.5 * m_size));
	int r = (dir.y > 0.9999f ? m_size - 1 : floor((dir.y + 1.0) * 0.5 * m_size));

	double *data = m_values->getFloat64Data();
	int idx = (m_size - r - 1) * m_size + c;
	for (int k = 0; k < 3; k++) {
		data[3 * idx + k] += (double)value[k] * dir.z * normFactor;
		m_totValue[k] += (double)value[k];
	}
	m_totValidParticles++;
}

void SphericalDistribution::scale(double scale) {
	double *data = m_values->getFloat64Data();
	for (int i = 0; i < m_size * m_size * SPECTRUM_SAMPLES; i++)
		*data++ *= scale;
}

void SphericalDistribution::saveExr(fs::path filename) {
	filename.replace_extension(".exr");
	ref<Bitmap> bitmap;
	bitmap = m_values->convert(Bitmap::ERGB, Bitmap::EFloat32);
	ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
	bitmap->write(Bitmap::EOpenEXR, stream);
}

void SphericalDistribution::load(Stream *stream) {
	m_size = stream->readInt();
	m_totValue[0] = stream->readDouble();
	m_totValue[1] = stream->readDouble();
	m_totValue[2] = stream->readDouble();
	m_totValidParticles = stream->readInt();
	stream->readDoubleArray(m_values->getFloat64Data(), m_size * m_size * SPECTRUM_SAMPLES);
}

void SphericalDistribution::save(Stream *stream) const {
	stream->writeInt(m_size);
	stream->writeDouble(m_totValue[0]);
	stream->writeDouble(m_totValue[1]);
	stream->writeDouble(m_totValue[2]);
	stream->writeInt(m_totValidParticles);
	stream->writeDoubleArray(m_values->getFloat64Data(), m_size * m_size * SPECTRUM_SAMPLES);
}

std::string SphericalDistribution::toString() const {
	std::ostringstream oss;
	oss << "SphericalDistribution[" << endl
		<< "  size = " << m_size << endl
		<< "]";
	return oss.str();
}

MultiLobeDistribution::MultiLobeDistribution(int numLobes, int size) : m_numLobes(numLobes) {
	m_lobes.resize(numLobes);
	for (int i = 0; i < numLobes; i++) {
		m_lobes[i] = new SphericalDistribution(size);
	}
}

SphericalDistribution* MultiLobeDistribution::getLobe(int lobeIdx) {
	return m_lobes[lobeIdx].get();
}

const SphericalDistribution* MultiLobeDistribution::getLobe(int lobeIdx) const {
	return m_lobes[lobeIdx].get();
}

void MultiLobeDistribution::clear() {
	for (int i = 0; i < m_numLobes; i++) {
		m_lobes[i]->clear();
	}
}

void MultiLobeDistribution::put(const MultiLobeDistribution *dist) {
	for (int i = 0; i < m_numLobes; i++) {
		m_lobes[i]->put(dist->getLobe(i));
	}
}

void MultiLobeDistribution::load(Stream *stream) {
	m_numLobes = stream->readInt();
	m_lobes.resize(m_numLobes);
	for (int i = 0; i < m_numLobes; i++) {
		m_lobes[i]->load(stream);
	}
}

void MultiLobeDistribution::save(Stream *stream) const {
	stream->writeInt(m_numLobes);
	for (int i = 0; i < m_numLobes; i++) {
		m_lobes[i]->save(stream);
	}
}

std::string MultiLobeDistribution::toString() const {
	std::ostringstream oss;
	oss << "MultiLobeDistribution[" << endl
		<< "  numLobes = " << m_numLobes << endl
		<< "]";
	return oss.str();
}

// WorkProcessor
BSDFRayTracer::BSDFRayTracer(const Vector &wi, int sqrtNumParticles, int size, 
		const AABB2 &aabb, int minDepth, int maxDepth, int shadowOption)
	: m_wi(wi), m_sqrtNumParticles(sqrtNumParticles), m_size(size), m_minDepth(minDepth), m_maxDepth(maxDepth), 
	m_aabb(aabb), m_shadowOption(shadowOption) { }

ref<WorkUnit> BSDFRayTracer::createWorkUnit() const {
	return new RangeWorkUnit();
}

ref<WorkResult> BSDFRayTracer::createWorkResult() const {
	return new MultiLobeDistribution(m_minDepth + 2, m_size);
}

ref<WorkProcessor> BSDFRayTracer::clone() const {
	return new BSDFRayTracer(m_wi, m_sqrtNumParticles, m_size, m_aabb, m_minDepth, m_maxDepth, m_shadowOption);
}

void BSDFRayTracer::prepare() {
	m_scene = static_cast<Scene *>(getResource("scene"));
	m_sampler = static_cast<Sampler *>(getResource("sampler"));
}

void BSDFRayTracer::process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop) {
	const RangeWorkUnit *range = static_cast<const RangeWorkUnit *>(workUnit);
	MultiLobeDistribution *res = static_cast<MultiLobeDistribution *>(workResult);

	double normFactor = (double)m_size * m_size / ((double)m_sqrtNumParticles * (double)m_sqrtNumParticles) * 0.25;
	//m_sampler->generate(Point2i(0));

	res->clear();
	
	for (size_t i = range->getRangeStart(); i <= range->getRangeEnd() && !stop; i++) {
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
		Spectrum throughput = sampleReflectance(ray, rRec, its);

		if (throughput[0] < -Epsilon)
			continue;
		//Assert(ray.d.z > -Epsilon);
		if (ray.d.z < 0)
			throughput = Spectrum(0.f);

		if (!throughput.isZero()) {
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
				res->getLobe(rRec.depth - 1)->put(ray.d, throughput, normFactor);
			else
				res->getLobe(m_minDepth)->put(ray.d, throughput, normFactor);

			res->getLobe(m_minDepth + 1)->put(ray.d, throughput, normFactor);
		}
	}
}

Point BSDFRayTracer::sampleRayOrigin(int idx, bool &success) {
	int r = idx / m_sqrtNumParticles;
	int c = idx % m_sqrtNumParticles;
	double x = m_aabb.min.x + (c + m_sampler->next1D()) / (double)m_sqrtNumParticles * (m_aabb.max.x - m_aabb.min.x);
	double y = m_aabb.min.y + (r + m_sampler->next1D()) / (double)m_sqrtNumParticles * (m_aabb.max.y - m_aabb.min.y);

	/*
	Point o(x, y, 1e2);
	Ray ray(o, Vector(0, 0, -1.0f), 0);
	Intersection its;
	m_scene->rayIntersect(ray, its);

	ray = Ray(its.p + m_wi * Epsilon, m_wi, 0);
	success = !m_scene->rayIntersect(ray);

	Point res = its.p + m_wi * 1e3;
	return res;
	*/

	Point o(x, y, 10.0f);
	o += m_wi * 1e3;
	Ray ray(o, -m_wi, 0);
	success = m_scene->rayIntersect(ray);
	
	//Intersection its;
	//success = m_scene->rayIntersect(ray, its);
	//Log(EInfo, "%d", m_aabb.contains(Point2(x, y)));
	//Log(EInfo, "(%.6f, %.6f), (%.6f, %.6f, %.6f)", x, y,
	//	its.p.x, its.p.y, its.p.z);

	return o;
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
			//if (phaseVal == 0)
			//	break;

			// Trace a ray
			ray = Ray(mRec.p, pRec.wo, ray.time);
			ray.mint = 0;
			scene->rayIntersect(ray, its);
		}
		else {
			if (rRec.medium)
				throughput *= mRec.transmittance / mRec.pdfFailure;

			if (!its.isValid() || !m_aabb.contains(Point2(its.p.x, its.p.y))) {
			//if (!its.isValid()) {
				//if (rRec.depth == 0)
				//	Log(EInfo, "%d, (%.6f, %.6f, %.6f), %d", rRec.depth, its.p.x, its.p.y, its.p.z,
				//		m_aabb.contains(Point2(its.p.x, its.p.y)));
				return throughput;
			}

			const BSDF *bsdf = its.getBSDF();
			BSDFSamplingRecord bRec(its, rRec.sampler);

			// Mark back-faced intersection as invalid
			//if (Frame::cosTheta(bRec.wi) <= 0) {
			//	return Spectrum(-1.0f);
			//}

			Spectrum bsdfVal = bsdf->sample(bRec, rRec.nextSample2D());
			throughput *= bsdfVal;
			//if (bsdfVal.isZero())
			//	break;

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
	m_maxDepth = stream->readInt();
	m_aabb = AABB2(stream);
}

void BSDFRayTracer::serialize(Stream *stream, InstanceManager *manager) const {
	for (int i = 0; i < 3; i++)
		stream->writeFloat(m_wi[i]);
	stream->writeInt(m_sqrtNumParticles);
	stream->writeInt(m_size);
	stream->writeInt(m_maxDepth);
	m_aabb.serialize(stream);
}

// ParallelProcess
BSDFSimulatorProcess::BSDFSimulatorProcess(const Vector &wi, int sqrtNumParticles,
		int size, const AABB2 &aabb, int minDepth, int maxDepth, int shadowOption)
		: m_wi(wi), m_sqrtNumParticles(sqrtNumParticles), m_size(size), m_minDepth(minDepth),
		m_maxDepth(maxDepth), m_aabb(aabb), m_shadowOption(shadowOption) {
	m_numParticles = m_sqrtNumParticles * m_sqrtNumParticles;
	m_start = 0;
	m_granularity = std::max((size_t)1, m_numParticles 
		/ (16 * Scheduler::getInstance()->getWorkerCount()));
	
	m_res = new MultiLobeDistribution(m_minDepth + 2, m_size);
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

	m_res->put(res);

	m_progress->update(++m_resultCount);
	UniqueLock lock(m_resultMutex);
	lock.unlock();
}

ref<WorkProcessor> BSDFSimulatorProcess::createWorkProcessor() const {
	return new BSDFRayTracer(m_wi, m_sqrtNumParticles, m_size, m_aabb, m_minDepth, m_maxDepth, m_shadowOption);
}

MTS_IMPLEMENT_CLASS(SphericalDistribution, false, WorkResult)
MTS_IMPLEMENT_CLASS(MultiLobeDistribution, false, WorkResult)
MTS_IMPLEMENT_CLASS(BSDFRayTracer, false, WorkProcessor)
MTS_IMPLEMENT_CLASS(BSDFSimulatorProcess, false, ParallelProcess)
MTS_NAMESPACE_END
