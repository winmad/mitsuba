#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/plugin.h>
#include "microfacet.h"

MTS_NAMESPACE_BEGIN

class RoughBSDFSingle: public BSDF {
public:
	RoughBSDFSingle(const Properties &props) : BSDF(props) {
		// avoid negative value
		m_offset = 1e4;

		ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

		Float uvscale = props.getFloat("uvscale", 1.0f);
		m_uvScale = Vector2(
			props.getFloat("uscale", uvscale),
			props.getFloat("vscale", uvscale)
		);

		m_moments0Filename = props.getString("moments0Filename", "");
		m_moments1Filename = props.getString("moments1Filename", "");

		if (m_moments0Filename == "") {
// 			m_moments0 = new ConstantSpectrumTexture(
// 				props.getSpectrum("moments0", Spectrum(m_offset)));
		}

		// alpha_x = sigma_x * sqrt(0.5)
		// alpha_y = sigma_y * sqrt(0.5)

		if (m_moments1Filename == "") {
// 			Spectrum defaultMoments1;
// 			defaultMoments1[0] = m_offset + 0.5;
// 			defaultMoments1[1] = m_offset + 0.5;
// 			defaultMoments1[2] = m_offset;
// 			m_moments1 = new ConstantSpectrumTexture(
// 				props.getSpectrum("moments1", defaultMoments1));
		}

		m_sampleVisibility = props.getBoolean("sampleVisibility", true);
	}

	RoughBSDFSingle(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		// avoid negative value
		m_offset = 1e4;

		m_sampleVisibility = stream->readBool();
		m_uvScale = Vector2(stream->readFloat(), stream->readFloat());
		m_moments0Filename = stream->readString();
		m_moments1Filename = stream->readString();

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeBool(m_sampleVisibility);
		stream->writeFloat(m_uvScale.x);
		stream->writeFloat(m_uvScale.y);
		stream->writeString(m_moments0Filename);
		stream->writeString(m_moments1Filename);
	}

	void configure() {
		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide | ESpatiallyVarying);

		m_usesRayDifferentials = false;

		m_samplers.resize(233);
		m_samplers[0] = static_cast<Sampler *>(PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));
		m_samplers[0]->configure();
		for (int i = 1; i < 233; i++) {
			m_samplers[i] = m_samplers[0]->clone();
		}

		Log(EInfo, "Start loading moments");
		// load moments
		m_moments0 = new Bitmap(fs::path(m_moments0Filename));
		m_moments1 = new Bitmap(fs::path(m_moments1Filename));
		m_size = m_moments0->getSize();

// 		ref<Bitmap> img = m_moments1;
// 		Spectrum spec = img->getPixel(Point2i(112, 111));
// 		Log(EInfo, "%.6f, %.6f, %.6f", spec[0], spec[1], spec[2]);
// 		/* 10156.716797 */
// 		exit(0);

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

// 		Log(EInfo, "wiLocal = (%.6f, %.6f, %.6f)", bRec.wi.x, bRec.wi.y, bRec.wi.z);
// 		Log(EInfo, "wiMacro = (%.6f, %.6f, %.6f)", wiMacro.x, wiMacro.y, wiMacro.z);

		Point2 uv = transformUV(bRec.its.uv);
		Point2i texP(math::floorToInt(uv.x * m_size.x), math::floorToInt(uv.y * m_size.y));
		Spectrum moments0 = m_moments0->getPixel(texP) - Spectrum(m_offset);
		Spectrum moments1 = m_moments1->getPixel(texP) - Spectrum(m_offset);
		Float sigmaX2 = moments1[0] - moments0[0] * moments0[0];
		Float sigmaY2 = moments1[1] - moments0[1] * moments0[1];
		Float cxy = moments1[2] - moments0[0] * moments0[1];
		Float sigmaX = std::sqrt(std::max(1e-8, sigmaX2));
		Float sigmaY = std::sqrt(std::max(1e-8, sigmaY2));
		Float rxy = cxy / (sigmaX * sigmaY);

		Normal mesoN(-moments0[0], -moments0[1], 1.0);
		mesoN = normalize(mesoN);

		Float cosWiMeso = dot(wiMacro, mesoN);
		Float cosWoMeso = dot(woMacro, mesoN);
		if (cosWiMeso <= 0 || cosWoMeso <= 0)
			return Spectrum(0.0f);

		ref<Sampler> sampler = m_samplers[Thread::getID() % 233];
		ref<BSDF> bsdf = m_bsdf;

		// no importance sampling... bad...

		// sample the slope
		Point2 sample2 = sampler->next2D();
		Point2 stdSample2 = warp::squareToStdNormal(sample2);		
		Vector2 slope;
		slope.x = stdSample2[0] * sigmaX + moments0[0];
		slope.y = (rxy * stdSample2[0] + std::sqrt(std::max(1e-8, 1.0 - rxy * rxy)) * stdSample2[1]) * sigmaY + moments0[1];

		Normal wm(-slope.x, -slope.y, 1.0);
// 		if (wm.length() < Epsilon) {
// 			Log(EInfo, "=================");
// 			Log(EInfo, "uv = (%.6f, %.6f)", bRec.its.uv.x * m_uvScale.x, bRec.its.uv.y * m_uvScale.y);
// 			Log(EInfo, "transformed_uv = (%.6f, %.6f)", uv.x, uv.y);
// 			Log(EInfo, "m0 = (%.6f, %.6f)", moments0[0], moments0[1]);
// 			Log(EInfo, "m1 = (%.6f, %.6f, %.6f)", moments1[0], moments1[1], moments1[2]);
// 			Log(EInfo, "cov matrix = (%.6f, %.6f; %.6f, %.6f)", sigmaX2, cxy, cxy, sigmaY2);
// 			Log(EInfo, "det = %.6f", sigmaX2 * sigmaY2 - cxy * cxy);
// 			Log(EInfo, "cov = (%.6f, %.6f, %.6f)", sigmaX, sigmaY, rxy);
// 			Log(EInfo, "%.6f, %.6f", sample2.x, sample2.y);
// 			Log(EInfo, "%.6f, %.6f", stdSample2[0], stdSample2[1]);
// 			Log(EInfo, "%.6f, %.6f", slope.x, slope.y);
// 		}
		wm = normalize(wm);

		Frame nFrame(wm);
		BSDFSamplingRecord bsdfRec(bRec.its, nFrame.toLocal(wiMacro), nFrame.toLocal(woMacro));
		Spectrum spec = bsdf->eval(bsdfRec);
		Spectrum res = spec * std::max(0.0, dot(wm, wiMacro)) / Frame::cosTheta(wm);
		res *= Frame::cosTheta(mesoN) / std::max(1e-4, dot(wiMacro, mesoN));
		
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

			res *= G2_G1;
		}

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
		if (m_bsdf->getClass()->getName() == "SmoothDiffuse")
			bRec.sampledType = EDiffuseReflection;
		else
			bRec.sampledType = EGlossyReflection;

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
		if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
			m_bsdf = static_cast<BSDF *>(child);
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
		oss << "RoughBSDFMulti[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  sampleVisibility " << m_sampleVisibility << "," << endl
			<< "  moments0 = " << indent(m_moments0->toString()) << "," << endl
			<< "  moments1 = " << indent(m_moments1->toString()) << "," << endl
			<< "]";
		return oss.str();
	}

	inline Point2 transformUV(const Point2 &_uv) const {
		Point2 uv(_uv);
		uv.x *= m_uvScale.x;
		uv.y *= m_uvScale.y;
		uv.x = uv.x - math::floorToInt(uv.x);
		uv.y = uv.y - math::floorToInt(uv.y);
		return uv;
	}

	MTS_DECLARE_CLASS()
private:
	ref<Bitmap> m_moments0, m_moments1;
	Vector2i m_size;

	Float m_offset;
	bool m_sampleVisibility;
	ref<BSDF> m_bsdf;

	std::string m_moments0Filename, m_moments1Filename;
	Vector2 m_uvScale;

	ref_vector<Sampler> m_samplers;
};

MTS_IMPLEMENT_CLASS_S(RoughBSDFSingle, false, BSDF)
MTS_EXPORT_PLUGIN(RoughBSDFSingle, "LEADR: Rough general BRDF with anisotropic Gaussian/GGX distribution");
MTS_NAMESPACE_END