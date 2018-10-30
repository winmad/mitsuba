#include <mitsuba/core/fstream.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include "scaleMatrixGrad_proc.h"

MTS_NAMESPACE_BEGIN

// WorkProcessor
ScaleMatrixGradTracer::ScaleMatrixGradTracer(const Vector &wi, int sqrtNumParticles, int wiResolution, 
		int woResolution, const AABB2 &aabb, int maxDepth, int shadowOption, 
		Bitmap *weightMap, Float distGiTexelScale)
	: m_wi(wi), m_sqrtNumParticles(sqrtNumParticles), m_wiResolution(wiResolution), m_woResolution(woResolution),
	m_maxDepth(maxDepth), m_aabb(aabb), m_shadowOption(shadowOption), 
	m_weightMap(weightMap), m_distGiTexelScale(distGiTexelScale) { }

ref<WorkUnit> ScaleMatrixGradTracer::createWorkUnit() const {
	return new RangeWorkUnit();
}

ref<WorkResult> ScaleMatrixGradTracer::createWorkResult() const {
	return new VectorBlock((2 * m_wiResolution * m_wiResolution * 2 * m_woResolution * m_woResolution + 1) * 3);
}

ref<WorkProcessor> ScaleMatrixGradTracer::clone() const {
	return new ScaleMatrixGradTracer(m_wi, m_sqrtNumParticles, m_wiResolution, m_woResolution, 
		m_aabb, m_maxDepth, m_shadowOption, m_weightMap, m_distGiTexelScale);
}

void ScaleMatrixGradTracer::prepare() {
	Scene *scene = static_cast<Scene *>(getResource("scene"));
	m_scene = new Scene(scene);
	m_sampler = static_cast<Sampler *>(getResource("sampler"));
	m_scene->setSampler(m_sampler);
}

void ScaleMatrixGradTracer::process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop) {
	const RangeWorkUnit *range = static_cast<const RangeWorkUnit *>(workUnit);
	VectorBlock *res = static_cast<VectorBlock *>(workResult);

	double normFactor = (double)m_woResolution * m_woResolution 
		/ (2.0 * M_PI * (double)m_sqrtNumParticles * m_sqrtNumParticles);
	
	//m_sampler->generate(Point2i(0));
	res->clear();

	m_distGiRange.x = (m_aabb.max.x - m_aabb.min.x) * m_distGiTexelScale;
	m_distGiRange.y = (m_aabb.max.y - m_aabb.min.y) * m_distGiTexelScale;

	std::vector<int> smIndices;
	std::vector<Float> smSumDiffs;
	
	// put throughput to res
	std::vector<int> fIndices(3);
	std::vector<Float> fValues(3);
	for (int k = 0; k < 3; k++) {
		fIndices[k] = (2 * m_wiResolution * m_wiResolution * 2 * m_woResolution * m_woResolution) * 3 + k;
	}
	
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
		smIndices.clear();
		smIndices.reserve(m_maxDepth * 60);
		smSumDiffs.clear();
		smSumDiffs.reserve(m_maxDepth * 60);
		Spectrum throughput = sampleReflectance(ray, rRec, its, smIndices, smSumDiffs);

// 		if (throughput[0] < -Epsilon) {
// 			//Log(EInfo, "invalid throughput");
// 			continue;
// 		}

		//Assert(ray.d.z > -Epsilon);
		if (ray.d.z < 0)
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

		//Log(EInfo, "Finish tracing one path.");

		// put grad contributions
		if (!throughput.isZero()) {
			Point2 wiTex = warp::uniformHemisphereToSquareConcentric(m_wi);
			Point2 woTex = warp::uniformHemisphereToSquareConcentric(ray.d);

			int woNumCells = m_woResolution - 1;
			int c1 = math::clamp(math::floorToInt(wiTex.x * m_wiResolution), 0, m_wiResolution - 1);
			int r1 = math::clamp(math::floorToInt(wiTex.y * m_wiResolution), 0, m_wiResolution - 1);
			int c2 = math::clamp(math::floorToInt(woTex.x * woNumCells), 0, woNumCells - 1);
			int r2 = math::clamp(math::floorToInt(woTex.y * woNumCells), 0, woNumCells - 1);

			Float u = woTex.x * woNumCells - c2;
			Float v = woTex.y * woNumCells - r2;

			Spectrum tmp = throughput * weight * normFactor;

			for (int dr = 0; dr < 2; dr++) {
				double wv = std::abs(1.0 - dr - v);
				for (int dc = 0; dc < 2; dc++) {
					double wKernel = wv * std::abs(1.0 - dc - u);
					
					Spectrum wGrad = m_weightMap->getPixel(Point2i(c2 + dc, r2 + dr));
					//Log(EInfo, "wGrad = %.6f, %.6f, %.6f", wGrad[0], wGrad[1], wGrad[2]);

					std::vector<Float> smGrads(smSumDiffs);
					for (int j = 0; j < smIndices.size(); j++) {
						int k = smIndices[j] % 3;
						smGrads[j] *= tmp[k] * wKernel * wGrad[k];
					}
					res->put(smIndices, smGrads, 0);
				}
			}

			for (int k = 0; k < 3; k++) {
				fValues[k] = throughput[k] * weight / 
					((double)m_sqrtNumParticles * (double)m_sqrtNumParticles);
			}
			res->put(fIndices, fValues, 0);

			//Log(EInfo, "after put res");
		}

		res->addWeight(weight);
	}
	//Log(EInfo, "process done %d: (%d, %d)", Thread::getID(), range->getRangeStart(), range->getRangeEnd());
}

Point ScaleMatrixGradTracer::sampleRayOrigin(int idx, double &weight) {
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
	}

	return o;
}

Spectrum ScaleMatrixGradTracer::sampleReflectance(RayDifferential &ray, RadianceQueryRecord &rRec, Intersection &getIts,
	std::vector<int> &smIndices, std::vector<Float> &smGrads) {
	const Scene *scene = rRec.scene;
	Intersection &its = rRec.its;
	MediumSamplingRecord mRec;

	rRec.rayIntersect(ray);
	Spectrum throughput(1.0f);

	while (rRec.depth <= m_maxDepth) {
		getIts = its;
		if (throughput.isZero())
			break;

		if (!its.isValid() || !m_aabb.contains(Point2(its.p.x, its.p.y)) ||
			rRec.depth == m_maxDepth) {
// 			if (rRec.depth == 0)
// 				Log(EInfo, "%d, (%.6f, %.6f, %.6f), %d", rRec.depth, its.p.x, its.p.y, its.p.z,
// 					m_aabb.contains(Point2(its.p.x, its.p.y)));
			return throughput;
		}

		const BSDF *bsdf = its.getBSDF();
		BSDFSamplingRecord bRec(its, rRec.sampler);

		// Mark back-faced intersection as invalid
		//if (Frame::cosTheta(bRec.wi) <= 0) {
		//	return Spectrum(-1.0f);
		//}
	
// 		if (rRec.depth > 0) {
// 			Vector wiMacro = its.baseFrame.toLocal(-ray.d);
// 			if (wiMacro.z > 0) {
// 				Log(EInfo, "Macro = %.6f, %.6f, %.6f", wiMacro.x, wiMacro.y, wiMacro.z);
// 				Log(EInfo, "Local = %.6f, %.6f, %.6f", its.wi.x, its.wi.y, its.wi.z);
// 			}
// 		}
			
		Spectrum bsdfVal = bsdf->sample(bRec, rRec.nextSample2D());
		throughput *= bsdfVal;
		if (bsdfVal.isZero()) {
			//Log(EInfo, "zero bsdf, wo = (%.6f, %.6f, %.6f)", bRec.wo.x, bRec.wo.y, bRec.wo.z);
			break;
		}

		for (int i = 0; i < 16; i++) {
			for (int c = 0; c < 3; c++) {
				smIndices.push_back(bRec.smIndices[i] * 3 + c);
				smGrads.push_back(bRec.smDiff[i][c]);
			}
		}

		const Vector wo = its.toWorld(bRec.wo);
		ray = Ray(its.p, wo, ray.time);

		scene->rayIntersect(ray, its);

		rRec.depth++;
	}

	return throughput;
}

ScaleMatrixGradTracer::ScaleMatrixGradTracer(Stream *stream, InstanceManager *manager) {
	for (int i = 0; i < 3; i++)
		m_wi[i] = stream->readFloat();
	m_sqrtNumParticles = stream->readInt();
	m_wiResolution = stream->readInt();
	m_woResolution = stream->readInt();
	m_maxDepth = stream->readInt();
	m_aabb = AABB2(stream);
	m_shadowOption = stream->readInt();
	// potential error: m_weightMap missing!
	m_distGiTexelScale = stream->readFloat();
}

void ScaleMatrixGradTracer::serialize(Stream *stream, InstanceManager *manager) const {
	for (int i = 0; i < 3; i++)
		stream->writeFloat(m_wi[i]);
	stream->writeInt(m_sqrtNumParticles);
	stream->writeInt(m_wiResolution);
	stream->writeInt(m_woResolution);
	stream->writeInt(m_maxDepth);
	m_aabb.serialize(stream);
	stream->writeInt(m_shadowOption);
	// potential error: m_weightMap missing!
	stream->writeFloat(m_distGiTexelScale);
}

// ParallelProcess
ScaleMatrixGradProcess::ScaleMatrixGradProcess(const Vector &wi, int sqrtNumParticles,
		int wiResolution, int woResolution, const AABB2 &aabb, int maxDepth, int shadowOption, 
		Bitmap *weightMap, Float distGiTexelScale)
		: m_wi(wi), m_sqrtNumParticles(sqrtNumParticles), m_wiResolution(wiResolution), m_woResolution(woResolution),
		m_maxDepth(maxDepth), m_aabb(aabb), m_shadowOption(shadowOption), m_weightMap(weightMap),
		m_distGiTexelScale(distGiTexelScale) {
	m_numParticles = m_sqrtNumParticles * m_sqrtNumParticles;
	m_start = 0;
	m_granularity = std::max((size_t)1, m_numParticles 
		/ (4 * Scheduler::getInstance()->getWorkerCount()));
	
	m_res = new VectorBlock((2 * m_wiResolution * m_woResolution * 2 * m_wiResolution * m_woResolution + 1) * 3);
	m_res->clear();

	m_resultCount = 0;
	m_progress = new ProgressReporter("Rendering", (m_numParticles - 1) / m_granularity + 1, NULL);
	m_resultMutex = new Mutex();
}

ScaleMatrixGradProcess::~ScaleMatrixGradProcess() {
	if (m_progress)
		delete m_progress;
}

void ScaleMatrixGradProcess::setGranularity(int granularity) {
	m_granularity = granularity;
}

ParallelProcess::EStatus ScaleMatrixGradProcess::generateWork(WorkUnit *unit, int worker) {
	RangeWorkUnit &range = *static_cast<RangeWorkUnit *>(unit);

	if (m_start >= m_numParticles)
		return EFailure;

	int ed = std::min(m_start + m_granularity - 1, m_numParticles - 1);
	range.setRange(m_start, ed);
	m_start = ed + 1;

	return ESuccess;
}

void ScaleMatrixGradProcess::processResult(const WorkResult *result, bool cancelled) {
	if (cancelled)
		return;

	const VectorBlock *res = static_cast<const VectorBlock *>(result);
	UniqueLock lock(m_resultMutex);

	m_res->put(res);
	m_progress->update(++m_resultCount);

	lock.unlock();
}

ref<WorkProcessor> ScaleMatrixGradProcess::createWorkProcessor() const {
	return new ScaleMatrixGradTracer(m_wi, m_sqrtNumParticles, m_wiResolution, m_woResolution, 
		m_aabb, m_maxDepth, m_shadowOption,
		m_weightMap, m_distGiTexelScale);
}

MTS_IMPLEMENT_CLASS_S(ScaleMatrixGradTracer, false, WorkProcessor)
MTS_IMPLEMENT_CLASS(ScaleMatrixGradProcess, false, ParallelProcess)
MTS_NAMESPACE_END
