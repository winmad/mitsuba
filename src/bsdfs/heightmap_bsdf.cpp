#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/vmf.h>
#include "microfacet.h"

MTS_NAMESPACE_BEGIN

// to be update
// default world frame with n=(0,0,1)

class HeightmapBSDF : public BSDF {
public:
	HeightmapBSDF(const Properties &props) : BSDF(props) {
		m_downsampleScale = props.getInteger("downsampleScale", 1);
		m_texSize = props.getInteger("texSize", 1024);
		m_tileX = props.getInteger("tileX", 1);
		m_tileY = props.getInteger("tileY", 1);
	}

	HeightmapBSDF(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		m_downsampleScale = stream->readInt();
		m_texSize = stream->readInt();
		m_tileX = stream->readInt();
		m_tileY = stream->readInt();

		m_hmap = static_cast<Shape *>(manager->getInstance(stream));
		m_bsdf = static_cast<BSDF *>(manager->getInstance(stream));

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeInt(m_downsampleScale);
		stream->writeInt(m_texSize);
		stream->writeInt(m_tileX);
		stream->writeInt(m_tileY);
		
		manager->serialize(stream, m_hmap.get());
		manager->serialize(stream, m_bsdf.get());
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

		m_uvBlockSize.x = (Float)m_downsampleScale / ((Float)m_texSize * m_tileX);
		m_uvBlockSize.y = (Float)m_downsampleScale / ((Float)m_texSize * m_tileY);

		BSDF::configure();
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EGlossyReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		Spectrum res(0.0f);
		Intersection its;

		ref<Sampler> sampler = m_samplers[Thread::getID() % 233];
		Vector wiWorld = bRec.its.toWorld(bRec.wi);
		Vector woWorld = bRec.its.toWorld(bRec.wo);

		Vector wiMacro = bRec.its.baseFrame.toLocal(wiWorld);
		Vector woMacro = bRec.its.baseFrame.toLocal(woWorld);

		Point2 uv;
		int bx = (int)std::floor(bRec.its.uv.x / m_uvBlockSize.x);
		int by = (int)std::floor(bRec.its.uv.y / m_uvBlockSize.y);
		uv.x = m_uvBlockSize.x * (bx + sampler->next1D());
		uv.y = m_uvBlockSize.y * (by + sampler->next1D());

		Normal norm;
		m_hmap->getPosAndNormal(uv, NULL, &norm);
		Frame nFrame(norm);
		BSDFSamplingRecord bsdfRec(bRec.its, nFrame.toLocal(wiMacro), nFrame.toLocal(woMacro));
		res = m_bsdf->eval(bsdfRec);
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
		if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
			m_bsdf = static_cast<BSDF *>(child);
		}
		else if (child->getClass()->derivesFrom(MTS_CLASS(Shape))) {
			Log(EInfo, "load heightmap");
			m_hmap = static_cast<Shape *>(child);
		}
		else {
			BSDF::addChild(name, child);
		}
	}

	Float getRoughness(const Intersection &its, int component) const {
		return m_bsdf->getRoughness(its, component);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HeightmapBSDF[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  baseBSDF = " << m_bsdf->toString() << "," << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	int m_downsampleScale;
	int m_texSize;
	int m_tileX;
	int m_tileY;

	Vector2 m_uvBlockSize;
	ref<Shape> m_hmap;
	ref<BSDF> m_bsdf;
	ref_vector<Sampler> m_samplers;
};

MTS_IMPLEMENT_CLASS_S(HeightmapBSDF, false, BSDF)
MTS_EXPORT_PLUGIN(HeightmapBSDF, "Base BSDF convolve tabulated NDF (induced by a heightmap)");
MTS_NAMESPACE_END