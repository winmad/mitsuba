#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/plugin.h>
#include "microfacet.h"

MTS_NAMESPACE_BEGIN

class AnisotropicRoughDiffuse : public BSDF {
public:
	AnisotropicRoughDiffuse(const Properties &props) : BSDF(props) {
		// avoid negative value
		m_offset = 1e4;

		ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

		m_reflectance = new ConstantSpectrumTexture(props.getSpectrum(
			props.hasProperty("reflectance") ? "reflectance"
			: "diffuseReflectance", Spectrum(0.5f)));

		m_moments0 = new ConstantSpectrumTexture(
			props.getSpectrum("moments0", Spectrum(m_offset)));
		
		// alpha_x = sigma_x * sqrt(0.5)
		// alpha_y = sigma_y * sqrt(0.5)
		Spectrum defaultMoments1;
		defaultMoments1[0] = m_offset + 0.5;
		defaultMoments1[1] = m_offset + 0.5;
		defaultMoments1[2] = m_offset;
		m_moments1 = new ConstantSpectrumTexture(
			props.getSpectrum("moments1", defaultMoments1));

		m_sampleVisibility = props.getBoolean("sampleVisibility", true);
	}

	AnisotropicRoughDiffuse(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		// avoid negative value
		m_offset = 1e4;

		m_type = (MicrofacetDistribution::EType) stream->readUInt();
		m_sampleVisibility = stream->readBool();
		m_moments0 = static_cast<Texture *>(manager->getInstance(stream));
		m_moments1 = static_cast<Texture *>(manager->getInstance(stream));
		m_reflectance = static_cast<Texture *>(manager->getInstance(stream));

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeUInt((uint32_t)m_type);
		stream->writeBool(m_sampleVisibility);
		manager->serialize(stream, m_moments0.get());
		manager->serialize(stream, m_moments1.get());
		manager->serialize(stream, m_reflectance.get());
	}

	void configure() {
		m_reflectance = ensureEnergyConservation(m_reflectance, "reflectance", 1.0f);

		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide
			| ((!m_reflectance->isConstant() || !m_moments0->isConstant() || !m_moments1->isConstant())
			? ESpatiallyVarying : 0));

		m_usesRayDifferentials = m_reflectance->usesRayDifferentials() ||
			m_moments0->usesRayDifferentials() || m_moments1->usesRayDifferentials();

		m_samplers.resize(233);
		m_samplers[0] = static_cast<Sampler *>(PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));
		m_samplers[0]->configure();
		for (int i = 1; i < 233; i++) {
			m_samplers[i] = m_samplers[0]->clone();
		}

		BSDF::configure();
	}

	inline Float approxLambda(const Vector &w, const Spectrum &moments0,
		Float sigmaX2, Float sigmaY2, Float cxy) const {
		if (Frame::sinTheta(w) < Epsilon)
			return 0;

		Float cotTheta = Frame::cosTheta(w) / Frame::sinTheta(w);
		Float muPhi = Frame::cosPhi(w) * moments0[0] + Frame::sinPhi(w) * moments0[1];
		Float sigma2Phi = Frame::cosPhi2(w) * sigmaX2 + Frame::sinPhi2(w) * sigmaY2 +
			2.0 * Frame::cosPhi(w) * Frame::sinPhi(w) * cxy;
		Float v = (cotTheta - muPhi) / (std::sqrt(2.0 * sigma2Phi));

// 		Log(EInfo, "w = (%.6f, %.6f, %.6f)", w.x, w.y, w.z);
// 		Log(EInfo, "%.6f, %.6f, %.6f", cotTheta, muPhi, sigma2Phi);
// 		Log(EInfo, "v = %.6f", v);

		if (v < 0)
			return 1e8f;

		if (v < 1.6)
			return (1.0 - 1.259 * v + 0.396 * v * v) / (3.535 * v + 2.181 * v * v);
		else
			return 0;
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EGlossyReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		Vector wiWorld = bRec.its.toWorld(bRec.wi);
		Vector woWorld = bRec.its.toWorld(bRec.wo);

		Vector wiMacro = bRec.its.baseFrame.toLocal(wiWorld);
		Vector woMacro = bRec.its.baseFrame.toLocal(woWorld);
	
		Spectrum moments0 = m_moments0->eval(bRec.its) - Spectrum(m_offset);
		Spectrum moments1 = m_moments1->eval(bRec.its) - Spectrum(m_offset);
		Float sigmaX2 = moments1[0] - moments0[0] * moments0[0];
		Float sigmaY2 = moments1[1] - moments0[1] * moments0[1];
		Float cxy = moments1[2] - moments0[0] * moments0[1];
		Float sigmaX = std::sqrt(sigmaX2);
		Float sigmaY = std::sqrt(sigmaY2);
		Float rxy = cxy / (sigmaX * sigmaY);

		Normal mesoN(-moments0[0], -moments0[1], 1.0);
		mesoN = normalize(mesoN);

		if (dot(wiMacro, mesoN) <= 0)
			return Spectrum(0.0f);
		Spectrum res = INV_PI * Frame::cosTheta(mesoN) /
			std::max(0.0, dot(wiMacro, mesoN)) *
			m_reflectance->eval(bRec.its);
		if (res.isZero())
			return Spectrum(0.0f);

		// sample the slope
		Float r = 0.0;
		int sampleCnt = 1;

		ref<Sampler> sampler = m_samplers[Thread::getID() % 233];

		for (int i = 0; i < sampleCnt; i++) {
			Point2 sample2 = sampler->next2D();
			Point2 stdSample2 = warp::squareToStdNormal(sample2);		
			Vector2 slope;
			slope.x = stdSample2[0] * sigmaX + moments0[0];
			slope.y = (rxy * stdSample2[0] + std::sqrt(1.0 - rxy * rxy) * stdSample2[1]) * sigmaY + moments0[1];

			Normal wm(-slope.x, -slope.y, 1.0);
// 			if (wm.length() < Epsilon) {
// 				Log(EInfo, "=================");
// 				Log(EInfo, "uv = (%.6f, %.6f)", (bRec.its.uv.x - 0.5) * 4, (bRec.its.uv.y - 0.5) * 4);
// 				Log(EInfo, "m0 = (%.6f, %.6f)", moments0[0], moments0[1]);
// 				Log(EInfo, "m1 = (%.6f, %.6f, %.6f)", moments1[0], moments1[1], moments1[2]);
// 				Log(EInfo, "cov matrix = (%.6f, %.6f; %.6f, %.6f)", sigmaX2, cxy, cxy, sigmaY2);
// 				Log(EInfo, "det = %.6f", sigmaX2 * sigmaY2 - cxy * cxy);
// 				Log(EInfo, "cov = (%.6f, %.6f, %.6f)", sigmaX, sigmaY, rxy);
// 				Log(EInfo, "%.6f, %.6f", sample2.x, sample2.y);
// 				Log(EInfo, "%.6f, %.6f", stdSample2[0], stdSample2[1]);
// 				Log(EInfo, "%.6f, %.6f", slope.x, slope.y);
// 			}
			wm = normalize(wm);

			Float tmpR = std::max(0.0, dot(wm, wiMacro)) * std::max(0.0, dot(wm, woMacro)) /
				Frame::cosTheta(wm);

			// height-correlated Smith masking-shadowing function
			if (m_sampleVisibility) {
				//Float G1_wi = distr.smithG1(bRec.wi, wm);
				//Float G1_wo = distr.smithG1(bRec.wo, wm);
				//Float G2_G1 = G1_wo / (G1_wi + G1_wo - G1_wi*G1_wo);
				Float G2_G1 = 0.0;
				if (dot(wm, wiMacro) > Epsilon && dot(wm, woMacro) > Epsilon) {
					Float lambda_wi = approxLambda(wiMacro, moments0, sigmaX2, sigmaY2, cxy);
					Float lambda_wo = approxLambda(woMacro, moments0, sigmaX2, sigmaY2, cxy);
					G2_G1 = 1.0 / (1.0 + lambda_wi + lambda_wo);
				}

// 				Log(EInfo, "cosWi = %.6f, cosWo = %.6f", dot(wm, bRec.wi), dot(wm, bRec.wo));
// 				Log(EInfo, "G = %.6f", G2_G1);

				tmpR *= G2_G1;
			}

			r += tmpR;
		}

		res *= r / (Float)sampleCnt;
		return res;
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

		bRec.wo = warp::squareToCosineHemisphere(sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;

		pdf = warp::squareToCosineHemispherePdf(bRec.wo);
		return eval(bRec, ESolidAngle) / pdf;
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sp) const {
		if (!(bRec.typeMask & EGlossyReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);
		Float pdf;
		return sample(bRec, pdf, sp);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "moments0")
				m_moments0 = static_cast<Texture *>(child);
			else if (name == "moments1")
				m_moments1 = static_cast<Texture *>(child);
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
		return 1.0;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "AnisotropicRoughDiffuse[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  distribution = " << MicrofacetDistribution::distributionName(m_type) << "," << endl
			<< "  sampleVisibility " << m_sampleVisibility << "," << endl
			<< "  moments0 = " << indent(m_moments0->toString()) << "," << endl
			<< "  moments1 = " << indent(m_moments1->toString()) << "," << endl
			<< "  reflectance = " << indent(m_reflectance->toString()) << "," << endl
			<< "]";
		return oss.str();
	}

	//Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	MicrofacetDistribution::EType m_type;
	ref<Texture> m_reflectance;
	ref<Texture> m_moments0, m_moments1;
	Float m_offset;
	bool m_sampleVisibility;
	ref_vector<Sampler> m_samplers;
};


MTS_IMPLEMENT_CLASS_S(AnisotropicRoughDiffuse, false, BSDF)
MTS_EXPORT_PLUGIN(AnisotropicRoughDiffuse, "Rough diffuse BRDF with anisotropic Gaussian distribution");
MTS_NAMESPACE_END