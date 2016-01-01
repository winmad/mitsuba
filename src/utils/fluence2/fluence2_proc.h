#ifndef __FLUENCE2_PROC_H
#define __FLUENCE2_PROC_H

#include <mitsuba/core/sched.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/particleproc.h>
#include <mitsuba/render/range.h>
#include <mitsuba/render/volume.h>
#include "blockinfo.h"


MTS_NAMESPACE_BEGIN


class Fluence2WorkResult : public WorkResult
{
public:
    inline Fluence2WorkResult()
    {
        m_range = new RangeWorkUnit;
    }

    void load(Stream *stream);
    void save(Stream *stream) const;
    std::string toString() const;

    void clear();
    void putSurface(const Point &p, const Spectrum &w);
    void putVolume(const Point &p, const Vector &d, const Spectrum &w);

    inline const RangeWorkUnit* getRangeWorkUnit() const
    {
        return m_range.get();
    }
    inline void setRangeWorkUnit(const RangeWorkUnit *range)
    {
        m_range->set(range);
    }

    MTS_DECLARE_CLASS()

protected:
    virtual ~Fluence2WorkResult() {}

public:
    std::vector<Point> m_surfPos, m_volPos;
    std::vector<Vector> m_volDir;
    std::vector<Spectrum> m_surfWeight, m_volWeight;

protected:
    ref<RangeWorkUnit> m_range;
};


class Fluence2ParticleTracer : public ParticleTracer
{
public:
    enum EmitterMode
    {
        EMITTER_AREA = 0,
        EMITTER_POINT = 1
    };

    Fluence2ParticleTracer(int maxDepth, const VolumeDataSource *vol);
    Fluence2ParticleTracer(int maxDepth, const VolumeDataSource *vol, const AABB &aabb);
    Fluence2ParticleTracer(Stream *stream, InstanceManager *manager);

    void serialize(Stream *stream, InstanceManager *manager) const;

    ref<WorkProcessor> clone() const;
    ref<WorkResult> createWorkResult() const;

    void process(const WorkUnit *workUnit, WorkResult *workResult,
        const bool &stop);

    void handleSurfaceInteraction(int depth, int nullInteractions, bool caustic,
        const Intersection &its, const Medium *medium,
        const Spectrum &weight);
    
    void handleMediumInteraction(int depth, int nullInteractions,
        bool delta, const MediumSamplingRecord &mRec, const Medium *medium,
        const Vector &wi, const Spectrum &weight);

    MTS_DECLARE_CLASS()

protected:
    virtual ~Fluence2ParticleTracer() {}

    const VolumeDataSource *m_densityVol;
    EmitterMode m_mode;
    AABB m_emitterAABB;
    Vector m_emitterExtents;
    ref<Fluence2WorkResult> m_result;
};


class Flunence2ParticleProcess : public ParticleProcess
{
public:
    Flunence2ParticleProcess(size_t sampleCount, size_t granularity, int maxDepth,
        const VolumeDataSource *vol, const AABB& medimAABB, const Vector3u &surfaceReso,
        const Vector3u &volumeReso, Float emitterArea, const char *msg = NULL);
    Flunence2ParticleProcess(size_t sampleCount, size_t granularity, int maxDepth,
        const VolumeDataSource *vol, const AABB& medimAABB, const Vector3u &surfaceReso,
        const Vector3u &volumeReso, const AABB &emitterAABB, const char *msg = NULL);

    virtual ~Flunence2ParticleProcess() {}

    void processResult(const WorkResult *wr, bool cancelled);
    ref<WorkProcessor> createWorkProcessor() const;

    inline const BlockInfo& getOutput() const
    {
        return m_output;
    }

    MTS_DECLARE_CLASS()

protected:
    void init(Float emitterArea, const AABB &volumeAABB, const Vector3u &volumeReso);

    int m_maxDepth;
    const VolumeDataSource *m_densityVol;
    Vector3u m_volReso;

    Fluence2ParticleTracer::EmitterMode m_emitterMode;
    AABB m_emitterAABB;
    
    BlockInfo m_output;

    Float m_surfScale, m_volScale;
};


MTS_NAMESPACE_END

#endif
