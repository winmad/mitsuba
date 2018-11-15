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

// to be update
// default world frame with n=(0,0,1)

class TabulatedBTF6D : public BSDF {
public:
	TabulatedBTF6D(const Properties &props) : BSDF(props) {
		m_xyReso = Vector2i(
			props.getInteger("xReso", 1),
			props.getInteger("yReso", 1));
		m_wiReso = props.getInteger("wiReso", 1);
		m_woReso = props.getInteger("woReso", 1);

		Float uvscale = props.getFloat("uvscale", 1.0f);
		m_uvScale = Vector2(
			props.getFloat("uscale", uvscale),
			props.getFloat("vscale", uvscale));

		m_filename = props.getString("filename", "");
	}

	TabulatedBTF6D(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		NotImplementedError("constructor");

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		NotImplementedError("serialize");
	}

	void configure() {
		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide | ESpatiallyVarying);

		m_usesRayDifferentials = false;

		m_btf = new Bitmap(fs::path(m_filename));
		Log(EInfo, "%d, %d", m_btf->getSize().x, m_btf->getSize().y);

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

		Intersection its(bRec.its);
		its.hasUVPartials = false;
		Point2 uv = transformUV(bRec.its.uv);
		int rowIdx = math::floorToInt(uv.y * m_xyReso.y) * m_xyReso.x + 
			math::floorToInt(uv.x * m_xyReso.x);

		Point2 wiTex = warp::uniformHemisphereToSquareConcentric(wiMacro);
		int wiNumCells = m_wiReso - 1;
		int c1 = math::clamp(math::floorToInt(wiTex.x * wiNumCells), 0, wiNumCells - 1);
		int r1 = math::clamp(math::floorToInt(wiTex.y * wiNumCells), 0, wiNumCells - 1);

		Point2 woTex = warp::uniformHemisphereToSquareConcentric(woMacro);
		int woNumCells = m_woReso - 1;
		int c2 = math::clamp(math::floorToInt(woTex.x * woNumCells), 0, woNumCells - 1);
		int r2 = math::clamp(math::floorToInt(woTex.y * woNumCells), 0, woNumCells - 1);

		Spectrum res(0.f);
		for (int dr1 = 0; dr1 < 2; dr1++) {
			Float v1 = wiTex.y * wiNumCells - r1;
			Float wv1 = std::abs(1.0 - dr1 - v1);

			for (int dc1 = 0; dc1 < 2; dc1++) {
				Float u1 = wiTex.x * wiNumCells - c1;
				Float wu1 = std::abs(1.0 - dc1 - u1);

				int wiIdx = (r1 + dr1) * m_wiReso + (c1 + dc1);
				for (int dr2 = 0; dr2 < 2; dr2++) {
					Float v2 = woTex.y * woNumCells - r2;
					Float wv2 = std::abs(1.0 - dr2 - v2);

					for (int dc2 = 0; dc2 < 2; dc2++) {
						Float u2 = woTex.x * woNumCells - c2;
						Float wu2 = std::abs(1.0 - dc2 - u2);

						int woIdx = (r2 + dr2) * m_woReso + (c2 + dc2);
						int colIdx = wiIdx * (m_woReso * m_woReso) + woIdx;
						Spectrum tmpValue = m_btf->getPixel(Point2i(colIdx, rowIdx));
						Float weight = wv1 * wu1 * wv2 * wu2;
						res += tmpValue * weight;
					}
				}
			}
		}

		//Log(EInfo, "wi = (%.6f, %.6f, %.6f), wo = (%.6f, %.6f, %.6f)", 
		//	bRec.wi.x, bRec.wi.y, bRec.wi.z, bRec.wo.x, bRec.wo.y, bRec.wo.z);

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

	inline Point2 transformUV(const Point2 &_uv) const {
		Point2 uv(_uv);
		uv.x *= m_uvScale.x;
		uv.y *= m_uvScale.y;
		uv.x = uv.x - math::floorToInt(uv.x);
		uv.y = uv.y - math::floorToInt(uv.y);
		return uv;
	}

	Float getRoughness(const Intersection &its, int component) const {
		return 0.f;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "TabulatedBTF[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  xResolution = " << m_xyReso.x << "," << endl
			<< "  yResolution = " << m_xyReso.y << "," << endl
			<< "  woResolution = " << m_woReso << "," << endl
			<< "  uScale = " << m_uvScale.x << "," << endl
			<< "  vScale = " << m_uvScale.y << "," << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Vector2i m_xyReso;
	int m_wiReso;
	int m_woReso;
	Vector2 m_uvScale;
	std::string m_filename;

	ref<Bitmap> m_btf;
};

MTS_IMPLEMENT_CLASS_S(TabulatedBTF6D, false, BSDF)
MTS_EXPORT_PLUGIN(TabulatedBTF6D, "Tabulated BTF");
MTS_NAMESPACE_END