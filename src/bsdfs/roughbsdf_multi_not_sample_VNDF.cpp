#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/plugin.h>
#include "microfacet.h"
#include "microfacet_multi.h"

#include <Eigen/Dense>

MTS_NAMESPACE_BEGIN

Spectrum samplePhaseFunction(const Intersection &its, Vector &wo, BSDF *bsdf, Sampler *sampler) {
	BSDFSamplingRecord bsdfRec(its, sampler);
	Spectrum res = bsdf->sample(bsdfRec, sampler->next2D());
	wo = bsdfRec.wo;
	return res;
}

Vector2 sampleSlopeGaussian(const Spectrum &moments0, Float sigmaX, Float sigmaY, Float rxy, Sampler *sampler) {
	Point2 sample2 = sampler->next2D();
	Point2 stdSample2 = warp::squareToStdNormal(sample2);
	Vector2 slope;
	slope.x = stdSample2[0] * sigmaX + moments0[0];
	slope.y = (rxy * stdSample2[0] + std::sqrt(1.0 - rxy * rxy) * stdSample2[1]) * sigmaY + moments0[1];
	return slope;
}

Vector2 sampleSlopeGGX(const Spectrum &moments0, Float sigmaX2, Float sigmaY2, Float cxy, Sampler *sampler) {
	Eigen::Matrix2d m;
	m(0, 0) = sigmaX2; m(0, 1) = cxy;
	m(1, 0) = cxy; m(1, 1) = sigmaY2;
	Eigen::JacobiSVD<Eigen::Matrix2d> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);

	Point2 sample2 = sampler->next2D();
	Float r = std::sqrt(sample2.x / (1.0 - sample2.x));
	Float phi = 2.0 * M_PI * sample2.y;
	Vector2 slope;
	slope.x = std::sqrt(2.0 * svd.singularValues()(0)) * r * std::cos(phi);
	slope.y = std::sqrt(2.0 * svd.singularValues()(1)) * r * std::sin(phi);

	Eigen::Matrix2d U = svd.matrixU();
	Float cosAngle = U(0, 0);
	Float sinAngle = U(1, 0);
	Vector2 res;
	res.x = slope.x * cosAngle - slope.y * sinAngle + moments0[0];
	res.y = slope.x * sinAngle + slope.y * cosAngle + moments0[1];
	return res;
}

Spectrum evalMulti(const Intersection &its, const Vector &wi, const Vector &wo, const Spectrum &moments0,
		Float sigmaX2, Float sigmaY2, Float cxy, BSDF *bsdf, int scatteringOrderMax, Sampler *sampler) {
	RayInfo ray;
	ray.updateDirection(-wi, moments0, sigmaX2, sigmaY2, cxy);
	//ray.updateDirection(-wi, 0.5, 0.5);

	RayInfo rayShadowing;
	rayShadowing.updateDirection(wo, moments0, sigmaX2, sigmaY2, cxy);
	//rayShadowing.updateDirection(wo, 0.5, 0.5);

	Spectrum res(0.0);
	
	ray.updateHeight(1.0);
	Spectrum energy(1.0);

	Float sigmaX = sqrt(sigmaX2);
	Float sigmaY = sqrt(sigmaY2);
	Float rxy = cxy / (sigmaX * sigmaY);

	int currentScatteringOrder = 0;
	while (currentScatteringOrder < scatteringOrderMax) {
		Float U = sampler->next1D();
		ray.updateHeight(sampleHeight(ray, U));

		if (ray.h == std::numeric_limits<Float>::max())
			break;
		else
			currentScatteringOrder++;

		// sample slope
		//Vector2 slope = sampleSlopeGaussian(moments0, sigmaX, sigmaY, rxy, sampler);
		Vector2 slope = sampleSlopeGGX(moments0, sigmaX2, sigmaY2, cxy, sampler);
		Normal wm(-slope.x, -slope.y, 1.0);
		wm = normalize(wm);

//		printf("=================\n");
// 		printf("uv = (%.6f, %.6f)\n", (its.uv.x - 0.5) * 4, (its.uv.y - 0.5) * 4);
// 		printf("m0 = (%.6f, %.6f)\n", moments0[0], moments0[1]);
// 		printf("cov matrix = (%.6f, %.6f; %.6f, %.6f)\n", sigmaX2, cxy, cxy, sigmaY2);
// 		printf("det = %.6f\n", sigmaX2 * sigmaY2 - cxy * cxy);
// 		printf("%.6f, %.6f\n", sample2.x, sample2.y);
// 		printf("%.6f, %.6f\n", stdSample2[0], stdSample2[1]);
// 		printf("%.6f, %.6f\n", slope.x, slope.y);

		// next event estimation
		Frame nFrame(wm);
		Float cosTerms = std::max(0.0, dot(-ray.w, wm)) / wm.z *
			(ray.mesoN.z / (ray.cosTheta * ray.Lambda));
		
		rayShadowing.updateHeight(ray.h);
		Float shadowingWo = rayShadowing.G1;
		BSDFSamplingRecord bsdfRec(its, nFrame.toLocal(-ray.w), nFrame.toLocal(wo));
		Spectrum spec = bsdf->eval(bsdfRec);
		Spectrum tmp = energy * bsdf->eval(bsdfRec) * cosTerms * shadowingWo;
		if (std::isfinite(tmp[0]) && std::isfinite(tmp[1]) && std::isfinite(tmp[2]))
			res += tmp;

		//if (tmp[0] < 0) {
		if (std::isinf(tmp[0])) {
		printf("=================\n");
		printf("mesoN = (%.6f, %.6f, %.6f)\n", ray.mesoN.x, ray.mesoN.y, ray.mesoN.z);
		printf("wi = (%.6f, %.6f, %.6f)\n", ray.w.x, ray.w.y, ray.w.z);
		printf("wo = (%.6f, %.6f, %.6f)\n", wo.x, wo.y, wo.z);
		printf("wm = (%.6f, %.6f, %.6f)\n", wm.x, wm.y, wm.z);
		printf("cosTerms = %.6f\n", cosTerms);
		printf("%.6f %.6f %.6f\n", dot(-ray.w, wm), ray.cosTheta, ray.Lambda);
		printf("bsdf spec = %.6f, %.6f, %.6f\n", spec[0], spec[1], spec[2]);
		printf("shadowing_wo = %.6f\n", shadowingWo);
		printf("cov matrix = (%.6f, %.6f; %.6f, %.6f)\n", sigmaX2, cxy, cxy, sigmaY2);
		printf("shadowing C1, Lambda = %.6f, %.6f\n", rayShadowing.C1, rayShadowing.Lambda);
 		printf("contrib = %.6f, %.6f, %.6f\n", tmp[0], tmp[1], tmp[2]);
		}

		// next direction
		Vector nextWr;
		Intersection tmpIts(its);
		tmpIts.wi = nFrame.toLocal(-ray.w);
		spec = samplePhaseFunction(tmpIts, nextWr, bsdf, sampler);
		if (spec.isZero())
			break;
		nextWr = nFrame.toWorld(nextWr);
		ray.updateDirection(nextWr, moments0, sigmaX2, sigmaY2, cxy);
		//ray.updateDirection(nextWr, 0.5, 0.5);
		energy = energy * spec * cosTerms;
		
// 		printf("=================\n");
// 		printf("weight = (%.6f, %.6f, %.6f)\n", spec[0], spec[1], spec[2]);
// 		printf("cosTerms = %.6f\n", cosTerms);
// 		printf("nextWr = (%.6f, %.6f, %.6f)\n", nextWr.x, nextWr.y, nextWr.z);
		
		if (!std::isfinite(energy[0]) || !std::isfinite(energy[1]) || !std::isfinite(energy[2]))
			break;
		ray.updateHeight(ray.h);

		// if NaN (should not happen, just in case)
		if ((ray.h != ray.h) || (ray.w.x != ray.w.x)) 
			return Spectrum(0.0f);
	}

	return res;
}

class RoughBSDFMulti : public BSDF {
public:
	RoughBSDFMulti(const Properties &props) : BSDF(props) {
		// avoid negative value
		m_offset = 5.0f;

		ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

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

		m_scatteringOrderMax = props.getInteger("scatteringOrderMax", 10);
	}

	RoughBSDFMulti(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
			// avoid negative value
			m_offset = 5.0f;

			m_sampleVisibility = stream->readBool();
			m_moments0 = static_cast<Texture *>(manager->getInstance(stream));
			m_moments1 = static_cast<Texture *>(manager->getInstance(stream));
			m_scatteringOrderMax = stream->readInt();

			configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeBool(m_sampleVisibility);
		manager->serialize(stream, m_moments0.get());
		manager->serialize(stream, m_moments1.get());
		stream->writeInt(m_scatteringOrderMax);
	}

	void configure() {
		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide
			| ((!m_moments0->isConstant() || !m_moments1->isConstant())
			? ESpatiallyVarying : 0));

		m_usesRayDifferentials = false;

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

		Float cosWiMeso = dot(wiMacro, mesoN);
		Float cosWoMeso = dot(woMacro, mesoN);
		if (cosWiMeso <= 0 || cosWoMeso <= 0)
			return Spectrum(0.0f);

		ref<Sampler> sampler = m_samplers[Thread::getID() % 233];
		ref<BSDF> bsdf = m_bsdf;

		Spectrum res;
		if (cosWiMeso < cosWoMeso) {
			res = evalMulti(bRec.its, wiMacro, woMacro, moments0, sigmaX2, sigmaY2, cxy, bsdf, m_scatteringOrderMax, sampler);
		}
		else {
			res = evalMulti(bRec.its, woMacro, wiMacro, moments0, sigmaX2, sigmaY2, cxy, bsdf, m_scatteringOrderMax, sampler) /
				Frame::cosTheta(wiMacro) * Frame::cosTheta(woMacro);
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
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "moments0")
				m_moments0 = static_cast<Texture *>(child);
			else if (name == "moments1")
				m_moments1 = static_cast<Texture *>(child);
			else
				BSDF::addChild(name, child);
		}
		else if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
			if (name == "baseBSDF")
				m_bsdf = static_cast<BSDF *>(child);
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
		oss << "RoughBSDFMulti[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  sampleVisibility " << m_sampleVisibility << "," << endl
			<< "  moments0 = " << indent(m_moments0->toString()) << "," << endl
			<< "  moments1 = " << indent(m_moments1->toString()) << "," << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_moments0, m_moments1;
	Float m_offset;
	bool m_sampleVisibility;
	int m_scatteringOrderMax;
	ref<BSDF> m_bsdf;

	ref_vector<Sampler> m_samplers;
};

MTS_IMPLEMENT_CLASS_S(RoughBSDFMulti, false, BSDF)
MTS_EXPORT_PLUGIN(RoughBSDFMulti, "Rough general BRDF with anisotropic Gaussian/GGX distribution");
MTS_NAMESPACE_END
