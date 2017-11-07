#pragma once
#if !defined(__BSDF_SIMULATOR_PROC_H_)
#define __BSDF_SIMULATOR_PROC_H_

#include <mitsuba/core/sched.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/aabb.h>
#include <mitsuba/render/range.h>
#include <boost/filesystem/path.hpp>

MTS_NAMESPACE_BEGIN

class SphericalDistribution : public WorkResult {
public:
	SphericalDistribution(int size);
	void clear();
	void put(const SphericalDistribution *dist);
	void put(const Vector &dir, const Spectrum &value, double normFactor);
	void scale(double scale);
	void saveExr(fs::path filename);

	void load(Stream *stream);
	void save(Stream *stream) const;
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	virtual ~SphericalDistribution() {}
public:
	int m_size;
	Vector3d m_totValue;
	int m_totValidParticles;
	ref<Bitmap> m_values;
};

class BSDFRayTracer : public WorkProcessor {
public:
	BSDFRayTracer(const Vector &wi, int sqrtNumParticles, int size, int maxDepth, const AABB2 &aabb);
	Point sampleRayOrigin(int idx);
	Spectrum sampleReflectance(RayDifferential &ray, RadianceQueryRecord &rRec);

	ref<WorkUnit> createWorkUnit() const;
	ref<WorkResult> createWorkResult() const;
	ref<WorkProcessor> clone() const;
	void prepare();
	void process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop);

	BSDFRayTracer(Stream *stream, InstanceManager *manager);
	void serialize(Stream *stream, InstanceManager *manager) const;

	MTS_DECLARE_CLASS()
protected:
	virtual ~BSDFRayTracer() {}
public:
	ref<Scene> m_scene;
	ref<Sampler> m_sampler;
	Vector m_wi;
	int m_sqrtNumParticles;
	int m_size;
	int m_maxDepth;
	AABB2 m_aabb;
};

class BSDFSimulatorProcess : public ParallelProcess {
public:
	BSDFSimulatorProcess(const Vector &wi, int sqrtNumParticles, int size, int maxDepth, const AABB2 &aabb);
	void setGranularity(int granularity);

	EStatus generateWork(WorkUnit *unit, int worker);
	void processResult(const WorkResult *result, bool cancelled);
	ref<WorkProcessor> createWorkProcessor() const;

	MTS_DECLARE_CLASS()
protected:
	virtual ~BSDFSimulatorProcess();
public:
	ref<Mutex> m_resultMutex;
	ProgressReporter *m_progress;
	int m_start;
	int m_numParticles;
	int m_granularity;
	int m_resultCount;

	Vector m_wi;
	int m_sqrtNumParticles;
	int m_size;
	int m_maxDepth;
	AABB2 m_aabb;
	ref<SphericalDistribution> m_res;
};

MTS_NAMESPACE_END

#endif /* __BSDF_SIMULATOR_PROC_H_ */