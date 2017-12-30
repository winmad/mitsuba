/*
	This plugin implements the BRDF model

	"Implementing a Simple Anisotropic Rough Diffuse
	Material with Stochastic Evaluation"

	from Eric Heitz and Jonathan Dupuy.
*/

#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>
#include "microfacet.h"

MTS_NAMESPACE_BEGIN

class SimpleAnisotropicRoughDiffuse : public BSDF {
public:
	SimpleAnisotropicRoughDiffuse(const Properties &props) : BSDF(props) {
		ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

		m_reflectance = new ConstantSpectrumTexture(props.getSpectrum(
			props.hasProperty("reflectance") ? "reflectance"
			: "diffuseReflectance", Spectrum(0.5f)));

		MicrofacetDistribution distr(props);
		m_type = distr.getType();
		m_sampleVisible = true; // this must be set to true

		m_alphaU = new ConstantFloatTexture(distr.getAlphaU());
		if (distr.getAlphaU() == distr.getAlphaV())
			m_alphaV = m_alphaU;
		else
			m_alphaV = new ConstantFloatTexture(distr.getAlphaV());
	}

	SimpleAnisotropicRoughDiffuse(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		m_type = (MicrofacetDistribution::EType) stream->readUInt();
		m_sampleVisible = stream->readBool();
		m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaV = static_cast<Texture *>(manager->getInstance(stream));
		m_reflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_sampleVisible = true; // this must be set to true

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeUInt((uint32_t)m_type);
		stream->writeBool(m_sampleVisible);
		manager->serialize(stream, m_alphaU.get());
		manager->serialize(stream, m_alphaV.get());
		manager->serialize(stream, m_reflectance.get());
	}

	void configure() {

		m_reflectance = ensureEnergyConservation(m_reflectance, "reflectance", 1.0f);

		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide
			| ((!m_reflectance->isConstant() || !m_alphaU->isConstant() || !m_alphaV->isConstant())
			? ESpatiallyVarying : 0));

		m_usesRayDifferentials = m_reflectance->usesRayDifferentials() ||
			m_alphaU->usesRayDifferentials() || m_alphaV->usesRayDifferentials();

		BSDF::configure();
	}

	/// Helper function: creates 1D uniform random numbers
	inline Float rand1D() const {
		return std::min(0.9999, Float(rand()) / Float(RAND_MAX));
	}

	/// Helper function: creates 2D uniform random numbers
	inline Point2 rand2D() const {
		return Point2(rand1D(), rand1D());
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EGlossyReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		// sample the VNDF
		MicrofacetDistribution distr(
			m_type,
			m_alphaU->eval(bRec.its).average(),
			m_alphaV->eval(bRec.its).average(),
			m_sampleVisible
			);
		Point2 sample2 = rand2D();
		Normal wm = distr.sample(bRec.wi, sample2);

		// height-correlated Smith masking-shadowing function
		Float G1_wi = distr.smithG1(bRec.wi, wm);
		Float G1_wo = distr.smithG1(bRec.wo, wm);
		Float G2_G1 = G1_wo / (G1_wi + G1_wo - G1_wi*G1_wo);

		return INV_PI * std::max(0.0, dot(wm, bRec.wo)) * G2_G1 * m_reflectance->eval(bRec.its);
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EGlossyReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return 0.0f;

		return warp::squareToCosineHemispherePdf(bRec.wo);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (!(bRec.typeMask & EGlossyReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);

		// sample the VNDF
		MicrofacetDistribution distr(
			m_type,
			m_alphaU->eval(bRec.its).average(),
			m_alphaV->eval(bRec.its).average(),
			m_sampleVisible
			);
		Normal wm = distr.sample(bRec.wi, sample);

		// diffuse reflection based on the microfacet normal

		// build orthonormal basis (Frisvad's method)
		Vector wx, wy;
		if (wm.z < -0.9999999f)
		{
			wx = Vector(0.0f, -1.0f, 0.0f);
			wy = Vector(-1.0f, 0.0f, 0.0f);
		}
		else {
			const Float a = 1.0f / (1.0f + wm.z);
			const Float b = -wm.x*wm.y*a;
			wx = Vector(1.0f - wm.x*wm.x*a, b, -wm.x);
			wy = Vector(b, 1.0f - wm.y*wm.y*a, -wm.y);
		}

		// diffuse reflection
		Point2 sample2 = rand2D();
		bRec.wo = warp::squareToCosineHemisphere(sample2);
		bRec.wo = wx * bRec.wo.x + wy * bRec.wo.y + wm * bRec.wo.z;
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;

		// side check
		if (Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		// height-correlated Smith masking-shadowing function
		Float G1_wi = distr.smithG1(bRec.wi, wm);
		Float G1_wo = distr.smithG1(bRec.wo, wm);
		Float G2_G1 = G1_wo / (G1_wi + G1_wo - G1_wi*G1_wo);
		Float weight = G2_G1;

		pdf = warp::squareToCosineHemispherePdf(bRec.wo);

		return weight * m_reflectance->eval(bRec.its);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sp) const {
		if (!(bRec.typeMask & EGlossyReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);
		Float pdf;
		return sample(bRec, pdf, sp);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "alpha")
				m_alphaU = m_alphaV = static_cast<Texture *>(child);
			else if (name == "alphaU")
				m_alphaU = static_cast<Texture *>(child);
			else if (name == "alphaV")
				m_alphaV = static_cast<Texture *>(child);
			else if (name == "reflectance" || name == "diffuseReflectance")
				m_reflectance = static_cast<Texture *>(child);
			else
				BSDF::addChild(name, child);
		}
		else {
			BSDF::addChild(name, child);
		}
	}

	Float getRoughness(const Intersection &its, int component) const {
		return 0.5f * (m_alphaU->eval(its).average()
			+ m_alphaV->eval(its).average());
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "AnisotropicRoughDiffuse[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  distribution = " << MicrofacetDistribution::distributionName(m_type) << "," << endl
			<< "  sampleVisible = " << m_sampleVisible << "," << endl
			<< "  alphaU = " << indent(m_alphaU->toString()) << "," << endl
			<< "  alphaV = " << indent(m_alphaV->toString()) << "," << endl
			<< "  reflectance = " << indent(m_reflectance->toString()) << "," << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	MicrofacetDistribution::EType m_type;
	ref<Texture> m_reflectance;
	ref<Texture> m_alphaU, m_alphaV;
	bool m_sampleVisible;
};

/**
* GLSL port of the rough conductor shader. This version is much more
* approximate -- it only supports the Ashikhmin-Shirley distribution,
* does everything in RGB, and it uses the Schlick approximation to the
* Fresnel reflectance of conductors. When the roughness is lower than
* \alpha < 0.2, the shader clamps it to 0.2 so that it will still perform
* reasonably well in a VPL-based preview.
*/
class SimpleAnisotropicRoughDiffuseShader : public Shader {
public:
	SimpleAnisotropicRoughDiffuseShader(Renderer *renderer, const Texture *reflectance, const Texture *alpha)
		: Shader(renderer, EBSDFShader), m_reflectance(reflectance), m_alpha(alpha) {
		m_reflectanceShader = renderer->registerShaderForResource(m_reflectance.get());
		m_alphaShader = renderer->registerShaderForResource(m_alpha.get());
	}

	bool isComplete() const {
		return m_reflectanceShader.get() != NULL &&
			m_alphaShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_reflectance.get());
		renderer->unregisterShaderForResource(m_alpha.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_reflectanceShader.get());
		deps.push_back(m_alphaShader.get());
	}

	void generateCode(std::ostringstream &oss,
		const std::string &evalName,
		const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) <= 0.0 || cosTheta(wo) <= 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    float sigma = " << depNames[1] << "(uv)[0] * 0.70711;" << endl
			<< "    float sigma2 = sigma * sigma;" << endl
			<< "    float A = 1.0 - 0.5 * sigma2 / (sigma2 + 0.33);" << endl
			<< "    float B = 0.45 * sigma2 / (sigma2 + 0.09);" << endl
			<< "    float maxCos = max(0.0, cosPhi(wi)*cosPhi(wo)+sinPhi(wi)*sinPhi(wo));" << endl
			<< "    float sinAlpha, tanBeta;" << endl
			<< "    if (cosTheta(wi) > cosTheta(wo)) {" << endl
			<< "        sinAlpha = sinTheta(wo);" << endl
			<< "        tanBeta = sinTheta(wi) / cosTheta(wi);" << endl
			<< "    } else {" << endl
			<< "        sinAlpha = sinTheta(wi);" << endl
			<< "        tanBeta = sinTheta(wo) / cosTheta(wo);" << endl
			<< "    }" << endl
			<< "    float value = A + B * maxCos * sinAlpha * tanBeta;" << endl
			<< "    return " << depNames[0] << "(uv) * inv_pi * cosTheta(wo) * value;" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) <= 0.0 || cosTheta(wo) <= 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return " << depNames[0] << "(uv) * inv_pi * cosTheta(wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_reflectance;
	ref<const Texture> m_alpha;
	ref<Shader> m_reflectanceShader;
	ref<Shader> m_alphaShader;
};

Shader *SimpleAnisotropicRoughDiffuse::createShader(Renderer *renderer) const {
	return new SimpleAnisotropicRoughDiffuseShader(renderer,
		m_reflectance.get(), m_alphaU.get());
}

MTS_IMPLEMENT_CLASS(SimpleAnisotropicRoughDiffuseShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(SimpleAnisotropicRoughDiffuse, false, BSDF)
MTS_EXPORT_PLUGIN(SimpleAnisotropicRoughDiffuse, "Centered axis-aligned anisotropic rough diffuse BRDF");
MTS_NAMESPACE_END