#include "fluence2_proc.h"
#include <mitsuba/core/warp.h>

// When this flag is on, the VPLs are uniformly sampled within each voxel.
// Otherwise, the samples at which the density values < 0.1 are rejected.
//#define VPL_FAST_MODE


MTS_NAMESPACE_BEGIN


void Fluence2WorkResult::load(Stream *stream)
{
    size_t sz;

    sz = stream->readSize();
    m_surfPos.resize(sz);
    m_surfWeight.resize(sz);
    stream->readArray<Point>(static_cast<Point *>(&m_surfPos[0]), sz);
    stream->readArray<Spectrum>(static_cast<Spectrum *>(&m_surfWeight[0]), sz);

    sz = stream->readSize();
    m_volPos.resize(sz);
    m_volWeight.resize(sz);
    stream->readArray<Point>(static_cast<Point *>(&m_volPos[0]), sz);
    stream->readArray<Spectrum>(static_cast<Spectrum *>(&m_volWeight[0]), sz);
}


void Fluence2WorkResult::save(Stream *stream) const
{
    size_t sz;

    sz = m_surfPos.size();
    stream->writeSize(sz);
    stream->writeArray<Point>(static_cast<const Point *>(&m_surfPos[0]), sz);
    stream->writeArray<Spectrum>(static_cast<const Spectrum *>(&m_surfWeight[0]), sz);

    sz = m_volPos.size();
    stream->writeSize(sz);
    stream->writeArray<Point>(static_cast<const Point *>(&m_volPos[0]), sz);
    stream->writeArray<Spectrum>(static_cast<const Spectrum *>(&m_volWeight[0]), sz);
}


std::string Fluence2WorkResult::toString() const
{
    std::ostringstream oss;
    oss << "Fluence2WorkResult[\n"
        << "    " << m_surfPos.size() << " surface points\n"
        << "    " << m_volPos.size() << " volume points\n"
        << "]";
    return oss.str();
}


void Fluence2WorkResult::clear()
{
    m_surfPos.clear();
    m_surfWeight.clear();

    m_volPos.clear();
    m_volDir.clear();
    m_volWeight.clear();
}


void Fluence2WorkResult::putSurface(const Point &p, const Spectrum &w)
{
    m_surfPos.push_back(p);
    m_surfWeight.push_back(w);
}


void Fluence2WorkResult::putVolume(const Point &p, const Vector &d, const Spectrum &w)
{
    m_volPos.push_back(p);
    m_volWeight.push_back(w);
}


Fluence2ParticleTracer::Fluence2ParticleTracer(int maxDepth, const VolumeDataSource *vol)
    : ParticleTracer(maxDepth, 10000, false), m_densityVol(vol), m_mode(EMITTER_AREA)
{
    m_emitterAABB.reset();
    m_emitterExtents = Vector(0.0f);
}


Fluence2ParticleTracer::Fluence2ParticleTracer(int maxDepth, const VolumeDataSource *vol, const AABB &aabb)
    : ParticleTracer(maxDepth, 10000, false), m_densityVol(vol), m_mode(EMITTER_POINT)
    , m_emitterAABB(aabb)
{
    m_emitterExtents = m_emitterAABB.getExtents();
}


Fluence2ParticleTracer::Fluence2ParticleTracer(Stream *stream, InstanceManager *manager)
    : ParticleTracer(stream, manager)
{
    m_mode = static_cast<EmitterMode>(stream->readInt());
    if ( m_mode == EMITTER_POINT )
    {
        m_emitterAABB = AABB(stream);
        m_emitterExtents = m_emitterAABB.getExtents();
    }
    else
    {
        m_emitterAABB.reset();
        m_emitterExtents = Vector(0.0f);
    }
}


void Fluence2ParticleTracer::serialize(Stream *stream, InstanceManager *manager) const
{
    ParticleTracer::serialize(stream, manager);
    stream->writeInt(m_mode);
    if ( m_mode == EMITTER_POINT )
        m_emitterAABB.serialize(stream);
}


ref<WorkProcessor> Fluence2ParticleTracer::clone() const
{
    return m_mode == EMITTER_POINT ? new Fluence2ParticleTracer(m_maxDepth, m_densityVol, m_emitterAABB)
        : new Fluence2ParticleTracer(m_maxDepth, m_densityVol);
}


ref<WorkResult> Fluence2ParticleTracer::createWorkResult() const
{
    return new Fluence2WorkResult();
}


void Fluence2ParticleTracer::process(const WorkUnit *workUnit, WorkResult *workResult,
    const bool &stop)
{
    const RangeWorkUnit *range = static_cast<const RangeWorkUnit *>(workUnit);
    m_result = static_cast<Fluence2WorkResult *>(workResult);
    m_result->clear();
    m_result->setRangeWorkUnit(range);

    /**
     * The following implementation is modified from ParticleTracer::process()
     */

    MediumSamplingRecord mRec;
    Intersection its;
    ref<Sensor> sensor    = m_scene->getSensor();
    bool needsTimeSample  = sensor->needsTimeSample();
    PositionSamplingRecord pRec(sensor->getShutterOpen()
        + 0.5f * sensor->getShutterOpenTime());
    Ray ray;

    m_sampler->generate(Point2i(0));

    for (size_t index = range->getRangeStart(); index <= range->getRangeEnd() && !stop; ++index) {
        m_sampler->setSampleIndex(index);

        /* Sample an emission */
        if (needsTimeSample)
            pRec.time = sensor->sampleTime(m_sampler->next1D());

        const Medium *medium;

        Spectrum power;
        Ray ray;

        if ( m_mode == EMITTER_AREA )
        {
            const Emitter *emitter = NULL;
            if (m_emissionEvents) {
                /* Sample the position and direction component separately to
                    generate emission events */
                power = m_scene->sampleEmitterPosition(pRec, m_sampler->next2D());
                emitter = static_cast<const Emitter *>(pRec.object);
                medium = emitter->getMedium();

                /* Forward the sampling event to the attached handler */
                handleEmission(pRec, medium, power);

                DirectionSamplingRecord dRec;
                power *= emitter->sampleDirection(dRec, pRec,
                        emitter->needsDirectionSample() ? m_sampler->next2D() : Point2(0.5f));
                ray.setTime(pRec.time);
                ray.setOrigin(pRec.p);
                ray.setDirection(dRec.d);
            } else {
                /* Sample both components together, which is potentially
                   faster / uses a better sampling strategy */

                power = m_scene->sampleEmitterRay(ray, emitter,
                    m_sampler->next2D(), m_sampler->next2D(), pRec.time);
                medium = emitter->getMedium();
                handleNewParticle();
            }
        }
        else
        {
            Assert(!m_emissionEvents);
            Assert(m_scene->getMedia().size() == 1);

            Point p;
#ifdef VPL_FAST_MODE
            p = m_emitterAABB.min;
            p.x += m_emitterExtents.x*m_sampler->next1D();
            p.y += m_emitterExtents.y*m_sampler->next1D();
            p.z += m_emitterExtents.z*m_sampler->next1D();
#else
            bool done = false;
            for ( int i = 0; i < 500; ++i )
            {
                p = m_emitterAABB.min;
                p.x += m_emitterExtents.x*m_sampler->next1D();
                p.y += m_emitterExtents.y*m_sampler->next1D();
                p.z += m_emitterExtents.z*m_sampler->next1D();
                if ( m_densityVol->lookupFloat(p) > 0.1f )
                {
                    done = true;
                    break;
                }
            }
            if ( !done )
            {
                Log(EWarn, "Failed to locate a non-empty voxel in %s", m_emitterAABB.toString().c_str());
                continue;
            }
#endif

            ray.setTime(0.0f);
            ray.setOrigin(p);
            ray.setDirection(warp::squareToUniformSphere(m_sampler->next2D()));

            power = Spectrum(1.0f);
            medium = m_scene->getMedia()[0];
        }

        int depth = 1, nullInteractions = 0;
        bool delta = false;

        Spectrum throughput(1.0f); // unitless path throughput (used for russian roulette)
        while (!throughput.isZero() && (depth <= m_maxDepth || m_maxDepth < 0)) {
            m_scene->rayIntersectAll(ray, its);

            /* ==================================================================== */
            /*                 Radiative Transfer Equation sampling                 */
            /* ==================================================================== */
            if (medium && medium->sampleDistance(Ray(ray, 0, its.t), mRec, m_sampler)) {
                /* Sample the integral
                  \int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
                */

                throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

                mRec.hasExtraInfo = true;
                mRec.extra = Point(ray.d);

                /* Forward the medium scattering event to the attached handler */
                handleMediumInteraction(depth, nullInteractions,
                        delta, mRec, medium, -ray.d, throughput*power);

                PhaseFunctionSamplingRecord pRec(mRec, -ray.d, EImportance);

                throughput *= medium->getPhaseFunction()->sample(pRec, m_sampler);
                delta = false;

                ray = Ray(mRec.p, pRec.wo, ray.time);
                ray.mint = 0;
            } else if (its.t == std::numeric_limits<Float>::infinity()) {
                /* There is no surface in this direction */
                break;
            } else {
                /* Sample
                    tau(x, y) (Surface integral). This happens with probability mRec.pdfFailure
                    Account for this and multiply by the proper per-color-channel transmittance.
                */
                if (medium)
                {
                    Spectrum tmp = mRec.transmittance / mRec.pdfFailure;
                    //if ( std::abs(tmp.average() - 1.0f) > Epsilon )
                    //    Log(EWarn, "%s", tmp.toString().c_str());
                    throughput *= tmp;
                }

                const BSDF *bsdf = its.getBSDF();

                /* Forward the surface scattering event to the attached handler */
                handleSurfaceInteraction(depth, nullInteractions, delta, its, medium, throughput*power);

                BSDFSamplingRecord bRec(its, m_sampler, EImportance);
                Spectrum bsdfWeight = bsdf->sample(bRec, m_sampler->next2D());
                if (bsdfWeight.isZero())
                    break;

                /* Prevent light leaks due to the use of shading normals -- [Veach, p. 158] */
                Vector wi = -ray.d, wo = its.toWorld(bRec.wo);
                Float wiDotGeoN = dot(its.geoFrame.n, wi),
                      woDotGeoN = dot(its.geoFrame.n, wo);
                if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
                    woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                    break;

                /* Keep track of the weight, medium and relative
                   refractive index along the path */
                throughput *= bsdfWeight;
                if (its.isMediumTransition())
                    medium = its.getTargetMedium(woDotGeoN);

                if (bRec.sampledType & BSDF::ENull)
                    ++nullInteractions;
                else
                    delta = bRec.sampledType & BSDF::EDelta;

#if 0
                /* This is somewhat unfortunate: for accuracy, we'd really want the
                   correction factor below to match the path tracing interpretation
                   of a scene with shading normals. However, this factor can become
                   extremely large, which adds unacceptable variance to output
                   renderings.

                   So for now, it is disabled. The adjoint particle tracer and the
                   photon mapping variants still use this factor for the last
                   bounce -- just not for the intermediate ones, which introduces
                   a small (though in practice not noticeable) amount of error. This
                   is also what the implementation of SPPM by Toshiya Hachisuka does.

                   Ultimately, we'll need better adjoint BSDF sampling strategies
                   that incorporate these extra terms */

                /* Adjoint BSDF for shading normals -- [Veach, p. 155] */
                throughput *= std::abs(
                    (Frame::cosTheta(bRec.wi) * woDotGeoN)/
                    (Frame::cosTheta(bRec.wo) * wiDotGeoN));
#endif

                ray.setOrigin(its.p);
                ray.setDirection(wo);
                ray.mint = Epsilon;
            }

            if (depth++ >= m_rrDepth) {
                /* Russian roulette: try to keep path weights equal to one,
                   Stop with at least some probability to avoid
                   getting stuck (e.g. due to total internal reflection) */

                Float q = std::min(throughput.max(), (Float) 0.95f);
                if (m_sampler->next1D() >= q)
                    break;
                throughput /= q;
            }
        }
    }

    /**
     * End of ParticleTracer::process()
     */

    m_result = NULL;
}


void Fluence2ParticleTracer::handleSurfaceInteraction(int depth, int nullInteractions,
    bool caustic, const Intersection &its, const Medium *medium, const Spectrum &weight)
{
    if ( its.isMediumTransition() && its.wi.z < 0.0f )
    {
        m_result->putSurface(its.p, weight);
    }
}


void Fluence2ParticleTracer::handleMediumInteraction(int depth,
    int nullInteractions, bool delta, const MediumSamplingRecord &mRec,
    const Medium *medium, const Vector &wi, const Spectrum &weight)
{
    Assert(mRec.hasExtraInfo);
    m_result->putVolume(mRec.p, Vector(mRec.extra), weight/mRec.sigmaS);
}


Flunence2ParticleProcess::Flunence2ParticleProcess(size_t sampleCount, size_t granularity,
    int maxDepth, const VolumeDataSource *vol, const AABB& medimAABB, const Vector3u &surfaceReso,
    const Vector3u &volumeReso, Float emitterArea, const char *msg)
    : ParticleProcess(ParticleProcess::ETrace, sampleCount, granularity,
        msg ? msg : "Rendering", NULL)
    , m_maxDepth(maxDepth), m_densityVol(vol), m_volReso(volumeReso)
    , m_emitterMode(Fluence2ParticleTracer::EMITTER_AREA)
    , m_output(medimAABB, surfaceReso, volumeReso)
{
    init(emitterArea, medimAABB, volumeReso);
}


Flunence2ParticleProcess::Flunence2ParticleProcess(size_t sampleCount, size_t granularity,
    int maxDepth, const VolumeDataSource *vol, const AABB& medimAABB, const Vector3u &surfaceReso,
    const Vector3u &volumeReso, const AABB &emitterAABB, const char *msg)
    : ParticleProcess(ParticleProcess::ETrace, sampleCount, granularity,
        msg ? msg : "Rendering", NULL)
    , m_maxDepth(maxDepth), m_densityVol(vol), m_volReso(volumeReso)
    , m_emitterMode(Fluence2ParticleTracer::EMITTER_POINT), m_emitterAABB(emitterAABB)
    , m_output(medimAABB, surfaceReso, volumeReso)
{
    init(0.0f, medimAABB, volumeReso);
}


void Flunence2ParticleProcess::init(Float emitterArea, const AABB &volumeAABB, const Vector3u &volumeReso)
{
    switch ( m_emitterMode )
    {
    case Fluence2ParticleTracer::EMITTER_AREA:
        m_surfScale = 1.0f/(emitterArea*M_PI*static_cast<Float>(m_workCount));
        //m_volScale = static_cast<Float>(volumeReso.x*volumeReso.y*volumeReso.z)/
        //    (static_cast<Float>(m_workCount)*volumeAABB.getVolume());
        m_volScale = static_cast<Float>(volumeReso.x*volumeReso.y*volumeReso.z)/
            (emitterArea*M_PI*static_cast<Float>(m_workCount)*volumeAABB.getVolume());
        break;

    case Fluence2ParticleTracer::EMITTER_POINT:
        m_surfScale = 1.0f/static_cast<Float>(m_workCount);
        m_volScale = static_cast<Float>(volumeReso.x*volumeReso.y*volumeReso.z)/
            (static_cast<Float>(m_workCount)*volumeAABB.getVolume());
        break;

    default:
        m_surfScale = m_volScale = 1.0f;
    }
}


void Flunence2ParticleProcess::processResult(const WorkResult *wr, bool cancelled)
{
    if ( !cancelled )
    {
        LockGuard lock(m_resultMutex);
        const Fluence2WorkResult *res = static_cast<const Fluence2WorkResult *>(wr);
        const RangeWorkUnit *range = res->getRangeWorkUnit();
        size_t npoint;
    
        npoint = res->m_surfPos.size();
        for ( size_t i = 0; i < npoint; ++i )
            m_output.addSurfacePoint(res->m_surfPos[i], res->m_surfWeight[i]*m_surfScale);

        npoint = res->m_volPos.size();
        for ( size_t i = 0; i < npoint; ++i )
            m_output.addVolumePoint(res->m_volPos[i], res->m_volWeight[i]*m_volScale);

        increaseResultCount(range->getSize());
    }
}


ref<WorkProcessor> Flunence2ParticleProcess::createWorkProcessor() const
{
    switch ( m_emitterMode )
    {
    case Fluence2ParticleTracer::EMITTER_AREA:
        return new Fluence2ParticleTracer(m_maxDepth, m_densityVol);
    case Fluence2ParticleTracer::EMITTER_POINT:
        return new Fluence2ParticleTracer(m_maxDepth, m_densityVol, m_emitterAABB);
    default:
        Assert(false);
        return NULL;
    }
}


MTS_IMPLEMENT_CLASS(Fluence2WorkResult, false, WorkResult)
MTS_IMPLEMENT_CLASS(Fluence2ParticleTracer, false, ParticleTracer)
MTS_IMPLEMENT_CLASS(Flunence2ParticleProcess, false, ParticleProcess)
MTS_NAMESPACE_END
