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

class MultiLobeDistribution : public WorkResult {
public:
	MultiLobeDistribution(int numLobes, int size);
	SphericalDistribution *getLobe(int lobeIdx);
	const SphericalDistribution *getLobe(int lobeIdx) const;
	void clear();
	void put(const MultiLobeDistribution *dist);

	void load(Stream *stream);
	void save(Stream *stream) const;
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	virtual ~MultiLobeDistribution() {}
public:
	int m_numLobes;
	std::vector<ref<SphericalDistribution> > m_lobes;
};

class BSDFRayTracer : public WorkProcessor {
public:
	BSDFRayTracer(const Vector &wi, int sqrtNumParticles, int size, const AABB2 &aabb,
		int minDepth, int maxDepth, int shadowOption);
	Point sampleRayOrigin(int idx, bool &success);
	Spectrum sampleReflectance(RayDifferential &ray, RadianceQueryRecord &rRec, Intersection &getIts);

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
	int m_minDepth;
	int m_maxDepth;
	AABB2 m_aabb;
	int m_shadowOption;
};

class BSDFSimulatorProcess : public ParallelProcess {
public:
	BSDFSimulatorProcess(const Vector &wi, int sqrtNumParticles, int size, const AABB2 &aabb, 
		int minDepth, int maxDepth, int shadowOption);
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
	int m_minDepth;
	int m_maxDepth;
	AABB2 m_aabb;
	int m_shadowOption;
	ref<MultiLobeDistribution> m_res;
};

MTS_NAMESPACE_END

#endif /* __BSDF_SIMULATOR_PROC_H_ */