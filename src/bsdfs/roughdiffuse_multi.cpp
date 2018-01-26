#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/plugin.h>
#include "microfacet.h"
#include "microfacet_multi.h"

MTS_NAMESPACE_BEGIN

void buildOrthonormalBasis(vec3& omega_1, vec3& omega_2, const vec3& omega_3) {
	if (omega_3.z < -0.9999999f) {
		omega_1 = vec3(0.0f, -1.0f, 0.0f);
		omega_2 = vec3(-1.0f, 0.0f, 0.0f);
	} 
	else {
		const Float a = 1.0f / (1.0f + omega_3.z);
		const Float b = -omega_3.x * omega_3.y * a ;
		omega_1 = vec3(1.0f - omega_3.x * omega_3. x * a , b , -omega_3.x);
		omega_2 = vec3(b, 1.0f - omega_3.y * omega_3.y * a, -omega_3.y);
	}
}

vec3 samplePhaseFunction_diffuse(const vec3& wm, Sampler *sampler) {
	// sample diffuse reflection
	vec3 w1, w2;
	buildOrthonormalBasis(w1, w2, wm);

	vec3 v = warp::squareToCosineHemisphere(sampler->next2D());
	vec3 wo = v.x*w1 + v.y*w2 + v.z*wm;

	return wo;
}

Spectrum eval_diffuse(const vec3& wi, const vec3& wo, 
		const Float alpha_x, const Float alpha_y, const Spectrum& albedo, 
		const int scatteringOrderMax, Sampler *sampler) {
	if(wi.z <= 0 || wo.z <= 0)
		return Spectrum(0.0f);

	// init
	RayInfo ray;
	ray.updateDirection(-wi, alpha_x, alpha_y);	

	RayInfo ray_shadowing;
	ray_shadowing.updateDirection(wo, alpha_x, alpha_y);

	Spectrum res(0.0f);

	ray.updateHeight(1.0f);
	Spectrum energy(1.0f);

	// random walk
	int current_scatteringOrder = 0;	
	while(current_scatteringOrder < scatteringOrderMax) {
		// next height
		Float U = sampler->next1D();
		ray.updateHeight(sampleHeight(ray, U));		

		// leave the microsurface?
		if( ray.h == std::numeric_limits<Float>::max() )
			break;
		else
			current_scatteringOrder++;

		// sample VNDF
		vec3 wm = sampleVNDF(-ray.w, alpha_x, alpha_y, sampler);

		// next event estimation
		Spectrum phasefunction = albedo * std::max(0.0, dot(wm, wo)) / M_PI;
		if (current_scatteringOrder == 1) { 
			// closed masking and shadowing (we compute G2 / G1 because G1 is already in the phase function)			
			Float G2_G1 = (1.0f + (-ray.Lambda-1.0f)) / (1.0f + (-ray.Lambda-1.0f) + ray_shadowing.Lambda);
			Spectrum I = energy * phasefunction * G2_G1;
			if (std::isfinite(I[0]))
				res += I;
		}
		else {
			Spectrum phasefunction = albedo * std::max(0.0, dot(wm, wo)) / M_PI;
			ray_shadowing.updateHeight(ray.h);
			Float shadowing = ray_shadowing.G1;
			Spectrum I = energy * phasefunction * shadowing;
			if (std::isfinite(I[0]))
				res += I;
		}

		// next direction
		ray.updateDirection(samplePhaseFunction_diffuse(wm, sampler), alpha_x, alpha_y);
		energy = energy * albedo;
		ray.updateHeight(ray.h);

		// if NaN (should not happen, just in case)
		if((ray.h != ray.h) || (ray.w.x != ray.w.x)) 
			return Spectrum(0.0f);
	}

	return res;
}

vec3 sample_diffuse(const vec3& wi, const Float alpha_x, const Float alpha_y, 
		const Spectrum& albedo, const int scatteringOrderMax, Spectrum& energy,
		Sampler *sampler) {
	energy = Spectrum(1.0f);

	// init
	RayInfo ray;
	ray.updateDirection(-wi, alpha_x, alpha_y);
	ray.updateHeight(1.0f);

	// random walk
	int current_scatteringOrder = 0;	
	while(true) {
		// next height
		Float U = sampler->next1D();
		ray.updateHeight(sampleHeight(ray, U));		

		// leave the microsurface?
		if( ray.h == std::numeric_limits<Float>::max())
			break;
		else
			current_scatteringOrder++;

		// sample VNDF
		vec3 wm = sampleVNDF(-ray.w, alpha_x, alpha_y, sampler);

		// next direction
		ray.updateDirection(samplePhaseFunction_diffuse(wm, sampler), alpha_x, alpha_y);
		energy = energy * albedo;
		ray.updateHeight(ray.h);

		// if NaN (should not happen, just in case)
		if ((ray.h != ray.h) || (ray.w.x != ray.w.x)) {
			energy = Spectrum(0.0f);
			return vec3(0,0,1);
		}

		if (current_scatteringOrder > scatteringOrderMax) {
			energy = Spectrum(0.0f);
			return vec3(0,0,1);
		}
	}

	return ray.w;
}

class RoughDiffuseMulti : public BSDF {
public:
	RoughDiffuseMulti(const Properties &props) : BSDF(props) {
		ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

		/* For better compatibility with other models, support both
		   'reflectance' and 'diffuseReflectance' as parameter names */
		m_reflectance = new ConstantSpectrumTexture(props.getSpectrum(
			props.hasProperty("reflectance") ? "reflectance"
				: "diffuseReflectance", Spectrum(1.0f)));

		// roughness
		Float alphaU, alphaV;
		if (props.hasProperty("alpha")) {
			alphaU = alphaV = props.getFloat("alpha");
			if (props.hasProperty("alphaU") || props.hasProperty("alphaV"))
				SLog(EError, "Microfacet model: please specify either 'alpha' or 'alphaU'/'alphaV'.");
		} else if (props.hasProperty("alphaU") || props.hasProperty("alphaV")) {
			if (!props.hasProperty("alphaU") || !props.hasProperty("alphaV"))
				SLog(EError, "Microfacet model: both 'alphaU' and 'alphaV' must be specified.");
			if (props.hasProperty("alpha"))
				SLog(EError, "Microfacet model: please specify either 'alpha' or 'alphaU'/'alphaV'.");
			alphaU = props.getFloat("alphaU");
			alphaV = props.getFloat("alphaV");
		}
		m_alphaU = new ConstantFloatTexture(alphaU);
		if (alphaU == alphaV)
			m_alphaV = m_alphaU;
		else
			m_alphaV = new ConstantFloatTexture(alphaV);

		// scattering order
		m_scatteringOrderMax = props.getInteger("scatteringOrderMax", 10);
	}

	RoughDiffuseMulti(Stream *stream, InstanceManager *manager)
			: BSDF(stream, manager) {
		m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaV = static_cast<Texture *>(manager->getInstance(stream));
		m_reflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_scatteringOrderMax = stream->readInt();

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_alphaU.get());
		manager->serialize(stream, m_alphaV.get());
		manager->serialize(stream, m_reflectance.get());
		stream->writeInt(m_scatteringOrderMax);
	}

	void configure() {
		unsigned int extraFlags = 0;
		if (m_alphaU != m_alphaV)
			extraFlags |= EAnisotropic;

		if (!m_alphaU->isConstant() || !m_alphaV->isConstant() ||
			!m_reflectance->isConstant())
			extraFlags |= ESpatiallyVarying;

		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide | extraFlags);

		/* Verify the input parameters and fix them if necessary */
		m_reflectance = ensureEnergyConservation(
			m_reflectance, "reflectance", 1.0f);

		m_usesRayDifferentials =
			m_alphaU->usesRayDifferentials() ||
			m_alphaV->usesRayDifferentials() ||
			m_reflectance->usesRayDifferentials();

		m_samplers.resize(233);
		m_samplers[0] = static_cast<Sampler *>(PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));
		m_samplers[0]->configure();
		for (int i = 1; i < 233; i++) {
			m_samplers[i] = m_samplers[0]->clone();
		}

		BSDF::configure();
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		/* Stop if this component was not requested */
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);
		
		Spectrum albedo = m_reflectance->eval(bRec.its);
		vec3 wi(bRec.wi.x, bRec.wi.y, bRec.wi.z);
		vec3 wo(bRec.wo.x, bRec.wo.y, bRec.wo.z);

		const Float alpha_x = std::max(m_alphaU->eval(bRec.its).average(), (Float) 1e-4f);
		const Float alpha_y = std::max(m_alphaV->eval(bRec.its).average(), (Float) 1e-4f);

		ref<Sampler> sampler = m_samplers[Thread::getID() % 233];

		// start random walks from lower direction for eval
		Spectrum res = (wi.z < wo.z) ?
			eval_diffuse(wi, wo, alpha_x, alpha_y, albedo, m_scatteringOrderMax, sampler) :
			eval_diffuse(wo, wi, alpha_x, alpha_y, albedo, m_scatteringOrderMax, sampler) / 
				Frame::cosTheta(bRec.wi)*Frame::cosTheta(bRec.wo);

		return res;
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return 0.0f;
		
		return 1.0f / M_PI * Frame::cosTheta(bRec.wo);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		Float pdf;
		return this->sample(bRec, pdf, sample);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;	

		vec3 wi(bRec.wi.x, bRec.wi.y, bRec.wi.z);

		const Float alpha_x = std::max(m_alphaU->eval(bRec.its).average(), (Float) 1e-4f);
		const Float alpha_y = std::max(m_alphaV->eval(bRec.its).average(), (Float) 1e-4f);

		Spectrum albedo = m_reflectance->eval(bRec.its);
		Spectrum energy;
		ref<Sampler> sampler = m_samplers[Thread::getID() % 233];
		
		vec3 wo = sample_diffuse(wi, alpha_x, alpha_y, albedo, m_scatteringOrderMax, energy, sampler);
		bRec.wo.x = wo.x;
		bRec.wo.y = wo.y;
		bRec.wo.z = wo.z;

		pdf = this->pdf(bRec, ESolidAngle);
		
		return energy;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "alpha")
				m_alphaU = m_alphaV = static_cast<Texture *>(child);
			else if (name == "alphaU")
				m_alphaU = static_cast<Texture *>(child);
			else if (name == "alphaV")
				m_alphaV = static_cast<Texture *>(child);
			else if (name == "reflectance")
				m_reflectance = static_cast<Texture *>(child);
			else
				BSDF::addChild(name, child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	Float getRoughness(const Intersection &its, int component) const {
		return 0.5f * (m_alphaU->eval(its).average()
			+ m_alphaV->eval(its).average());
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "RoughConductor[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  alphaU = " << indent(m_alphaU->toString()) << "," << endl
			<< "  alphaV = " << indent(m_alphaV->toString()) << "," << endl
			<< "  specularReflectance = " << indent(m_reflectance->toString()) << "," << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_reflectance;
	ref<Texture> m_alphaU, m_alphaV;
	int m_scatteringOrderMax;

	ref_vector<Sampler> m_samplers;
};

MTS_IMPLEMENT_CLASS_S(RoughDiffuseMulti, false, BSDF)
MTS_EXPORT_PLUGIN(RoughDiffuseMulti, "Rough diffuse (multi-scattering) BRDF")
MTS_NAMESPACE_END
