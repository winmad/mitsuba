#pragma once
#if !defined(__LOBE_DERIVATIVE_PROC_H_)
#define __LOBE_DERIVATIVE_PROC_H_

#include <mitsuba/core/sched.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/aabb.h>
#include <mitsuba/render/range.h>
#include <mitsuba/render/spherical_distribution.h>
#include <boost/filesystem/path.hpp>

MTS_NAMESPACE_BEGIN

class BSDFDerivativeRayTracer : public WorkProcessor {
public:
	BSDFDerivativeRayTracer(int numVars, const Vector &wi,
		int sqrtNumParticles, int size, const AABB2 &aabb,
		int minDepth, int maxDepth, int shadowOption);
	
	Point sampleRayOrigin(int idx, bool &success);
	void sampleDerivative(RayDifferential &ray, RadianceQueryRecord &rRec, 
		Intersection &getIts, Float normFactor, MultiLobeDistribution *res);

	ref<WorkUnit> createWorkUnit() const;
	ref<WorkResult> createWorkResult() const;
	ref<WorkProcessor> clone() const;
	void prepare();
	void process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop);

	BSDFDerivativeRayTracer(Stream *stream, InstanceManager *manager);
	void serialize(Stream *stream, InstanceManager *manager) const;

	MTS_DECLARE_CLASS()
protected:
	virtual ~BSDFDerivativeRayTracer() {}
public:
	int m_numVars;
	ref<Scene> m_scene;
	ref<Sampler> m_sampler;
	Vector m_wi;
	int m_sqrtNumParticles;
	int m_size;
	int m_minDepth;
	int m_maxDepth;
	AABB2 m_aabb;
	int m_shadowOption;
};

class BSDFDerivativeProcess : public ParallelProcess {
public:
	BSDFDerivativeProcess(int numVars, const Vector &wi,
		int sqrtNumParticles, int size, const AABB2 &aabb,
		int minDepth, int maxDepth, int shadowOption);
	void setGranularity(int granularity);

	EStatus generateWork(WorkUnit *unit, int worker);
	void processResult(const WorkResult *result, bool cancelled);
	ref<WorkProcessor> createWorkProcessor() const;

	MTS_DECLARE_CLASS()
protected:
	virtual ~BSDFDerivativeProcess();
public:
	ref<Mutex> m_resultMutex;
	ProgressReporter *m_progress;
	int m_start;
	int m_numParticles;
	int m_granularity;
	int m_resultCount;

	int m_numVars;
	Vector m_wi;
	int m_sqrtNumParticles;
	int m_size;
	int m_minDepth;
	int m_maxDepth;
	AABB2 m_aabb;
	int m_shadowOption;
	ref<MultiLobeDistribution> m_res;
};

MTS_NAMESPACE_END

#endif
