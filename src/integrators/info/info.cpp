#include <mitsuba/render/scene.h>
#include <mitsuba/render/lightcutter.h>
#include <mitsuba/render/vpl.h>

MTS_NAMESPACE_BEGIN
        
class InfoIntegrator : public SamplingIntegrator {
public:
    MTS_DECLARE_CLASS()
    std::string m_info_name;
    size_t m_samples;

    InfoIntegrator(const Properties &props) : SamplingIntegrator(props) {
        m_info_name = props.getString("infoName");
        m_samples = props.getSize("samples");
    }
    /// Unserialize from a binary data stream
    InfoIntegrator(Stream *stream, InstanceManager *manager)
            : SamplingIntegrator(stream, manager) {
        m_info_name = stream->readString();
        m_samples = stream->readSize();
    }
    
    void configureSampler(const Scene *scene, Sampler *sampler) {
        SamplingIntegrator::configureSampler(scene, sampler);
        if (m_samples > 1)
                sampler->request2DArray(m_samples);
    }
    
    void serialize(Stream *stream, InstanceManager *manager) const {
            SamplingIntegrator::serialize(stream, manager);
            stream->writeString(m_info_name);
            stream->writeSize(m_samples);
    }
    
    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        RayDifferential ray(r);
        Spectrum Li(0.0f);
        if (!rRec.rayIntersect(ray)) {
            return Spectrum(0.0f);
        }
        
        if(m_info_name == "normal") {
            Vector output_normal = its.shFrame.n+Vector(1.0f);
            return Spectrum((Float*)&output_normal);
        }
        
        if(m_info_name == "distance") {
            return Spectrum(its.t);
        }
        
        if(m_info_name == "filter_size") {
            Point2 sample;
            Point2 *sampleArray;
            if (m_samples > 1) {
                sampleArray = rRec.sampler->next2DArray(m_samples);
            }
            else {
                sample = rRec.nextSample2D(); sampleArray = &sample;
            }
            
            const BSDF *bsdf = its.getBSDF(ray);

            DirectSamplingRecord dRec(its);

            Float size = 1e3;
            
            if (bsdf->getType() & BSDF::ESmooth) {
                for (size_t i = 0; i<m_samples; ++i) {
                    scene->sampleEmitterDirect(dRec, sampleArray[i]);
                    Intersection occ_its;
                    Ray ray(dRec.ref, dRec.d, Epsilon,
                            dRec.dist*(1-ShadowEpsilon), dRec.time);
                    bool hit = scene->rayIntersect(ray, occ_its);
                    if(hit && occ_its.t < its.t) {
                        size = std::min((its.t/(its.t-occ_its.t)-1), size);
                    }
                }
            }
            return Spectrum(size);
        }
        
        return Li;
    }
};
        
MTS_IMPLEMENT_CLASS_S(InfoIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(InfoIntegrator, "Info integrator");
MTS_NAMESPACE_END