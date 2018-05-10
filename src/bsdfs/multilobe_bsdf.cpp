#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/vmf.h>
#include "microfacet.h"

MTS_NAMESPACE_BEGIN

// ************
// to be update
// default world frame with n=(0,0,1)
// ************
// UPDATE: use baseFrame
// ************
// WORKING: consider macro-surface deformation

class MultiLobeBSDF : public BSDF {
public:
	MultiLobeBSDF(const Properties &props) : BSDF(props) {
		m_numLobes = props.getInteger("numLobes", 1);
		
		Float uvscale = props.getFloat("uvscale", 1.0f);
		m_uvScale = Vector2(
			props.getFloat("uscale", uvscale),
			props.getFloat("vscale", uvscale)
			);

		m_lobeFilenamePrefix = props.getString("prefix", "");

		m_useMacroDeform = props.getBoolean("useMacroDeform", false);
		m_useApproxShadowing = props.getBoolean("useApproxShadowing", false);
	}

	MultiLobeBSDF(Stream *stream, InstanceManager *manager) 
		: BSDF(stream, manager) {
		m_numLobes = stream->readInt();
		m_uvScale = Vector2(stream->readFloat(), stream->readFloat());
		m_bsdf = static_cast<BSDF *>(manager->getInstance(stream));
		m_lobeFilenamePrefix = stream->readString();
		m_useMacroDeform = stream->readBool();
		m_useApproxShadowing = stream->readBool();

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeInt(m_numLobes);
		stream->writeFloat(m_uvScale.x);
		stream->writeFloat(m_uvScale.y);
		manager->serialize(stream, m_bsdf.get());
		stream->writeString(m_lobeFilenamePrefix);
		stream->writeBool(m_useMacroDeform);
		stream->writeBool(m_useApproxShadowing);
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

		// load textures
		m_lobes.resize(m_numLobes);
		for (int i = 0; i < m_numLobes; i++) {
			Properties props("bitmap");
			props.setString("wrapMode", "repeat");
			props.setFloat("gamma", 1.0f);
			props.setFloat("uscale", 1.0f);
			props.setFloat("vscale", 1.0f);
			props.setString("filterType", "nearest");
			std::ostringstream oss;
			oss << m_lobeFilenamePrefix << i << ".exr";
			props.setString("filename", oss.str());
			Log(EInfo, "Load vMF texture %s", oss.str().c_str());
			m_lobes[i] = static_cast<Texture *>(PluginManager::getInstance()->
				createObject(MTS_CLASS(Texture), props));
			m_lobes[i]->configure();
		}
		Log(EInfo, "Load vMF lobes as textures");

		/*
		Point2 uv = transformUV(Point2(0.6, 0.7));
		Log(EInfo, "%.6f, %.6f, %.6f", uv.x, uv.y * 0.5, uv.y * 0.5 + 0.5);
		
		Intersection its;
		its.uv = Point2(uv.x, uv.y * 0.5);
		Spectrum spec = m_lobes[0]->eval(its, false);
		Log(EInfo, "(%.6f, %.6f, %.6f)", spec[0], spec[1], spec[2]);
		its.uv = Point2(uv.x, uv.y * 0.5 + 0.5);
		spec = m_lobes[0]->eval(its, false);
		Log(EInfo, "(%.6f, %.6f, %.6f)", spec[0], spec[1], spec[2]);
		*/
		BSDF::configure();
	}

	inline Float conv(Float w1, Float lambda1, const Vector &mu1,
			Float w2, Float lambda2, const Vector &mu2) const {
		Vector d = mu1 * lambda1 + mu2 * lambda2;
		Float len = d.length();
		Float res = w1 * w2 * 4.0 * M_PI * math::fastexp(-lambda1 - lambda2) 
			* std::sinh(len) / len;
		
// 		if (!std::isfinite(res)) {
// 			std::ostringstream oss;
// 			oss << "======================\n"
// 				<< "w1 = " << w1 << ", " << "lambda1 = " << lambda1 << ", "
// 				<< "mu1 = (" << mu1.x << ", " << mu1.y << ", " << mu1.z << ")\n"
// 				<< "w2 = " << w2 << ", " << "lambda2 = " << lambda2 << ", "
// 				<< "mu2 = (" << mu2.x << ", " << mu2.y << ", " << mu2.z << ")\n";
// 			std::cout << oss.str();
// 		}
		
		return res;
	}

	Float G2(const Vector &wi, const Vector &wo, const Vector &nMeso, 
			const std::vector<Spectrum> &lobesParam0, 
			const std::vector<Spectrum> &lobesParam1) const {
		Float wCos = 1.1767;
		Float lambdaCos = 2.1440;
		Float integralWi = 0.0;
		Float integralWo = 0.0;
		for (int i = 0; i < m_numLobes; i++) {
			Float alpha = lobesParam0[i][0];
			if (alpha < 1e-4)
				continue;
			Float kappa = lobesParam0[i][1];
			Vector mu(2.0 * lobesParam1[i][0] - 1.0, 2.0 * lobesParam1[i][1] - 1.0, 2.0 * lobesParam1[i][2] - 1.0);
			Float wNDF = alpha * kappa / (2 * M_PI * (1 - math::fastexp(-2 * kappa)));

			integralWi += conv(wCos, lambdaCos, wi, wNDF, kappa, mu);
			integralWo += conv(wCos, lambdaCos, wo, wNDF, kappa, mu);
		}

		Float a = 1.0 / std::abs(nMeso.z);
		Float G1Wi = std::max(0.0, dot(wi, nMeso)) * a / integralWi;
		Float G1Wo = std::max(0.0, dot(wo, nMeso)) * a / integralWo;
		Float tmp = G1Wi * G1Wo;

// 		std::ostringstream oss;
// 		oss << "======================\n"
// 			<< "area_wi = " << std::max(0.0, dot(wi, nMeso)) * a << "\n"
// 			<< "area_wo = " << std::max(0.0, dot(wo, nMeso)) * a << "\n"
// 			<< "meso = (" << nMeso.x << ", " << nMeso.y << ", " << nMeso.z << ")\n"
// 			<< "integral_wi = " << integralWi << "\n"
// 			<< "integral_wo = " << integralWo << "\n";
// 		std::cout << oss.str();

		return tmp / (G1Wi + G1Wo + tmp);
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EGlossyReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		Spectrum res(0.0f);
		Intersection its(bRec.its);
		its.hasUVPartials = false;
		Point2 uv = transformUV(bRec.its.uv);
		
		ref<Sampler> sampler = m_samplers[Thread::getID() % 233];
		
		Vector wiWorld = bRec.its.toWorld(bRec.wi);
		Vector woWorld = bRec.its.toWorld(bRec.wo);
		
		Vector wiMacro = bRec.its.baseFrame.toLocal(wiWorld);
		Vector woMacro = bRec.its.baseFrame.toLocal(woWorld);

		// hack to avoid normal inconsistency
		//Vector wiWorld = bRec.wi;
		//Vector woWorld = bRec.wo;

		//std::ostringstream oss;

		Vector nMeso(0.0);
		std::vector<Spectrum> lobesParam0(m_numLobes);
		std::vector<Spectrum> lobesParam1(m_numLobes);

		for (int i = 0; i < m_numLobes; i++) {
			its.uv.x = uv.x;
			its.uv.y = uv.y * 0.5;
			Spectrum param0 = m_lobes[i]->eval(its, false);
			its.uv.x = uv.x;
			its.uv.y = uv.y * 0.5 + 0.5;
			Spectrum param1 = m_lobes[i]->eval(its, false);

			if (m_useApproxShadowing) {
				lobesParam0[i] = param0;
				lobesParam1[i] = param1;
			}

			Float alpha = param0[0];
			if (alpha < 1e-8)
				continue;

			Float kappa = param0[1];
			VonMisesFisherDistr vmf(kappa);
			Vector mu(2.0 * param1[0] - 1.0, 2.0 * param1[1] - 1.0, 2.0 * param1[2] - 1.0);
			nMeso += alpha * mu;

			Frame lobeFrame(mu);
			Vector norm = lobeFrame.toWorld(vmf.sample(Point2(sampler->next1D(), sampler->next1D())));
			Frame nFrame(norm);
			
			BSDFSamplingRecord bsdfRec(bRec.its, nFrame.toLocal(wiMacro), nFrame.toLocal(woMacro));
			Spectrum spec = m_bsdf->eval(bsdfRec);

// 			if (mu.length() < 1e-8 || wiWorld.length() < 1e-8 || woWorld.length() < 1e-8 ||
// 				norm.length() < 1e-8 || bsdfRec.wi.length() < 1e-8 || bsdfRec.wo.length() < 1e-8)
// 			oss << "===================" << std::endl
// 				<< "uv = " << uv.x << ", " << uv.y << std::endl
// 				<< "lobe " << i << std::endl
// 				<< "alpha = " << alpha << ", "
// 				<< "kappa = " << kappa << std::endl
// 				<< "mu = (" << mu.x << ", " << mu.y << ", " << mu.z << ")" << std::endl
// 				<< "len(mu) = " << mu.length() << std::endl
// 				<< "wi_world = (" << wiWorld.x << ", " << wiWorld.y << ", " << wiWorld.z << ")" << std::endl
// 				<< "wo_world = (" << woWorld.x << ", " << woWorld.y << ", " << woWorld.z << ")" << std::endl
// 				<< "normal = (" << norm.x << ", " << norm.y << ", " << norm.z << ")" << std::endl
// 				<< "wi_local = (" << bsdfRec.wi.x << ", " << bsdfRec.wi.y << ", " << bsdfRec.wi.z << ")" << std::endl
// 				<< "wo_local = (" << bsdfRec.wo.x << ", " << bsdfRec.wo.y << ", " << bsdfRec.wo.z << ")" << std::endl
// 				<< "value = (" << spec[0] << ", " << spec[1] << ", " << spec[2] << ")" << std::endl;

			res += alpha * spec;
		}

		if (m_useApproxShadowing) {
			Float len = nMeso.length();
			if (len < 1e-4)
				return Spectrum(0.0);
			nMeso /= len;

			Float shadowTerm = G2(wiMacro, woMacro, nMeso, lobesParam0, lobesParam1);
			if (std::isfinite(shadowTerm))
				res *= shadowTerm;
		}

		//std::cout << oss.str();

		return res;
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EGlossyReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return 0.0f;
		
		// hack...
		return warp::squareToCosineHemispherePdf(bRec.wo);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (!(bRec.typeMask & EGlossyReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);
		
		
		// naive sampling
		/*
		bRec.wo = warp::squareToCosineHemisphere(sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		pdf = warp::squareToCosineHemispherePdf(bRec.wo);
		return eval(bRec, ESolidAngle) / pdf;
		*/
		
		Spectrum res(0.0f);

		Intersection its(bRec.its);
		its.hasUVPartials = false;
		Point2 uv = transformUV(bRec.its.uv);

		Vector nMeso(0.0);
		std::vector<Spectrum> lobesParam0(m_numLobes);
		std::vector<Spectrum> lobesParam1(m_numLobes);

		std::vector<Float> cdf(0);
		std::vector<int> lobeIndices(0);
		for (int i = 0; i < m_numLobes; i++) {
			its.uv.x = uv.x;
			its.uv.y = uv.y * 0.5;
			Spectrum param0 = m_lobes[i]->eval(its, false);

			if (m_useApproxShadowing) {
				lobesParam0[i] = param0;
				lobesParam1[i] = m_lobes[i]->eval(its, false);
			}

			if (param0[0] < 1e-8)
				continue;
			
			lobeIndices.push_back(i);
			if (cdf.empty())
				cdf.push_back(param0[0]);
			else
				cdf.push_back(cdf.back() + param0[0]);
		}
		if (cdf.size() == 0) {
			Log(EError, "no lobes?!");
		}
		Assert(cdf.size() > 0);

		Float normFactor = cdf.back();
		for (int i = 0; i < cdf.size(); i++)
			cdf[i] /= normFactor;

		ref<Sampler> sampler = m_samplers[Thread::getID() % 233];
		int tmp = (int)(std::lower_bound(cdf.begin(), cdf.end(), sampler->next1D()) - cdf.begin());
		Assert(tmp < cdf.size() && tmp >= 0);
		int lobeIdx = lobeIndices[tmp];

		// sample a single vMF lobe
		Float alpha = cdf[lobeIdx];
		if (lobeIdx > 0)
			alpha -= cdf[lobeIdx - 1];

		its.uv.x = uv.x;
		its.uv.y = uv.y * 0.5;
		Spectrum param0 = m_lobes[lobeIdx]->eval(its, false);
		its.uv.x = uv.x;
		its.uv.y = uv.y * 0.5 + 0.5;
		Spectrum param1 = m_lobes[lobeIdx]->eval(its, false);

		Float kappa = param0[1];
		VonMisesFisherDistr vmf(kappa);
		Vector mu(2.0 * param1[0] - 1.0, 2.0 * param1[1] - 1.0, 2.0 * param1[2] - 1.0);

		Vector wiWorld = bRec.its.toWorld(bRec.wi);
		Vector wiMacro = bRec.its.baseFrame.toLocal(wiWorld);

		// hack to avoid normal inconsistency
		//Vector wiWorld = bRec.wi;

		Frame lobeFrame(mu);
		Vector norm = lobeFrame.toWorld(vmf.sample(Point2(sampler->next1D(), sampler->next1D())));
		Frame nFrame(norm);
		
		its = bRec.its;
		its.wi = nFrame.toLocal(wiMacro);
		BSDFSamplingRecord bsdfRec(its, sampler);

		res = m_bsdf->sample(bsdfRec, Point2(sampler->next1D(), sampler->next1D()));

		Vector woMacro = nFrame.toWorld(bsdfRec.wo);
		Vector woWorld = bRec.its.baseFrame.toWorld(woMacro);

		bRec.wo = bRec.its.toLocal(woWorld);
		bRec.sampledComponent = bsdfRec.sampledComponent;
		// diffuse or glossy
		bRec.sampledType = bsdfRec.sampledType;

		// hardcode: base brdf has only one lobe
		bRec.used[0] = true;
		bRec.dAlbedo[0] = bsdfRec.dAlbedo[0];
		bRec.dRoughness[0] = bsdfRec.dRoughness[0];
		bRec.used[1] = false;

		// hack..
		pdf = warp::squareToCosineHemispherePdf(bRec.wo);

		// side check
		if (Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		if (m_useApproxShadowing) {
			Float len = nMeso.length();
			if (len < 1e-4)
				return Spectrum(0.0);
			nMeso /= len;

			Float shadowTerm = G2(wiMacro, woMacro, nMeso, lobesParam0, lobesParam1);
			if (std::isfinite(shadowTerm))
				res *= shadowTerm;
		}
		
// 		std::ostringstream oss;
// 		oss << "===================" << std::endl
// 			<< "lobe " << lobeIdx << std::endl
// 			<< "alpha = " << alpha << ", "
// 			<< "kappa = " << kappa << std::endl
// 			<< "mu = (" << mu.x << ", " << mu.y << ", " << mu.z << ")" << std::endl;
// 		oss << "wi_world = (" << wiWorld.x << ", " << wiWorld.y << ", " << wiWorld.z << ")" << std::endl
// 			<< "wo_world = (" << woWorld.x << ", " << woWorld.y << ", " << woWorld.z << ")" << std::endl
// 			<< "normal = (" << norm.x << ", " << norm.y << ", " << norm.z << ")" << std::endl
// 			<< "wi_local_micro = (" << bsdfRec.wi.x << ", " << bsdfRec.wi.y << ", " << bsdfRec.wi.z << ")" << std::endl
// 			<< "wo_local_micro = (" << bsdfRec.wo.x << ", " << bsdfRec.wo.y << ", " << bsdfRec.wo.z << ")" << std::endl
// 			<< "value = (" << res[0] << ", " << res[1] << ", " << res[2] << ")" << std::endl
// 			<< "wi_macro = (" << wiMacro.x << ", " << wiMacro.y << ", " << wiMacro.z << ")" << std::endl
// 			<< "wo_macro = (" << woMacro.x << ", " << woMacro.y << ", " << woMacro.z << ")" << std::endl;
// 		std::cout << oss.str() << std::endl;

		return res;
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sp) const {
		if (!(bRec.typeMask & EGlossyReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);
		Float pdf;
		return sample(bRec, pdf, sp);
	}

	inline Point2 transformUV(const Point2 &_uv) const {
		Point2 uv(_uv);
		uv.x *= m_uvScale.x;
		uv.y *= m_uvScale.y;
		uv.x = uv.x - math::floorToInt(uv.x);
		uv.y = uv.y - math::floorToInt(uv.y);
		return uv;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
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
		return m_bsdf->getRoughness(its, component);
	}

	Spectrum getLobeAlbedo(const Intersection &its, int component) const {
		return m_bsdf->getLobeAlbedo(its, component);
	}

	Float getLobeRoughness(const Intersection &its, int component) const {
		return m_bsdf->getLobeRoughness(its, component);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MultiLobeBSDF[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  numLobes = " << m_numLobes  << "," << endl
			<< "  baseBSDF = " << m_bsdf->toString() << "," << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	int m_numLobes;
	std::string m_lobeFilenamePrefix;
	ref_vector<Texture> m_lobes;
	Vector2 m_uvScale;
	ref<BSDF> m_bsdf;
	ref_vector<Sampler> m_samplers;
	bool m_useMacroDeform;
	bool m_useApproxShadowing;
};

MTS_IMPLEMENT_CLASS_S(MultiLobeBSDF, false, BSDF)
MTS_EXPORT_PLUGIN(MultiLobeBSDF, "Base BSDF convolve multi-lobe NDF");
MTS_NAMESPACE_END