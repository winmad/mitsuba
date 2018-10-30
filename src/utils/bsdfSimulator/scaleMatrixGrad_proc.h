#pragma once
#if !defined(__SCALE_MATRIX_GRAD_PROC_H_)
#define __SCALE_MATRIX_GRAD_PROC_H_

#include <mitsuba/core/sched.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/aabb.h>
#include <mitsuba/render/range.h>
#include <mitsuba/render/vectorblock.h>
#include <boost/filesystem/path.hpp>

MTS_NAMESPACE_BEGIN

class ScaleMatrixGradTracer : public WorkProcessor {
public:
	ScaleMatrixGradTracer(const Vector &wi, int sqrtNumParticles, int wiResolution, int woResolution, const AABB2 &aabb,
		int maxDepth, int shadowOption, Bitmap *weightMap, Float distGiTexelScale);
	Point sampleRayOrigin(int idx, double &weight);
	Spectrum sampleReflectance(RayDifferential &ray, RadianceQueryRecord &rRec, Intersection &getIts,
		std::vector<int> &smIndices, std::vector<Float> &smGrads);

	ref<WorkUnit> createWorkUnit() const;
	ref<WorkResult> createWorkResult() const;
	ref<WorkProcessor> clone() const;
	void prepare();
	void process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop);

	ScaleMatrixGradTracer(Stream *stream, InstanceManager *manager);
	void serialize(Stream *stream, InstanceManager *manager) const;

	MTS_DECLARE_CLASS()
protected:
	virtual ~ScaleMatrixGradTracer() {}
public:
	ref<Scene> m_scene;
	ref<Sampler> m_sampler;
	Vector m_wi;
	int m_sqrtNumParticles;
	int m_wiResolution;
	int m_woResolution;
	int m_minDepth;
	int m_maxDepth;
	AABB2 m_aabb;
	int m_shadowOption;

	Bitmap *m_weightMap;
	Float m_distGiTexelScale;
	Vector2 m_distGiRange;
};

class ScaleMatrixGradProcess : public ParallelProcess {
public:
	ScaleMatrixGradProcess(const Vector &wi, int sqrtNumParticles, int wiResolution, int woResolution, const AABB2 &aabb, 
		int maxDepth, int shadowOption, Bitmap *weightMap, Float distGiTexelScale);
	void setGranularity(int granularity);

	EStatus generateWork(WorkUnit *unit, int worker);
	void processResult(const WorkResult *result, bool cancelled);
	ref<WorkProcessor> createWorkProcessor() const;

	MTS_DECLARE_CLASS()
protected:
	virtual ~ScaleMatrixGradProcess();
public:
	ref<Mutex> m_resultMutex;
	ProgressReporter *m_progress;
	int m_start;
	int m_numParticles;
	int m_granularity;
	int m_resultCount;

	Vector m_wi;
	int m_sqrtNumParticles;
	int m_wiResolution;
	int m_woResolution;
	int m_minDepth;
	int m_maxDepth;
	AABB2 m_aabb;
	int m_shadowOption;
	Bitmap *m_weightMap;
	Float m_distGiTexelScale;
	
	ref<VectorBlock> m_res;
};

MTS_NAMESPACE_END

#endif /* __BSDF_SIMULATOR_PROC_H_ */