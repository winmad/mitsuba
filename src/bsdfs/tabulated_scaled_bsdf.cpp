#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/mipmap.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

class TabulatedScaledBSDF : public BSDF {
public:
	TabulatedScaledBSDF(const Properties &props) : BSDF(props) {
		m_angularScaleFilename = props.getString("angularScaleFilename", "");
	}

	TabulatedScaledBSDF(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		m_angularScaleFilename = stream->readString();

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeString(m_angularScaleFilename);
	}

	void configure() {
		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide | ESpatiallyVarying);

		m_usesRayDifferentials = false;

		// load angular scales
		m_angularScales = new Bitmap(fs::path(m_angularScaleFilename));
		m_lobeSize = m_angularScales->getSize();
		m_wiResolution = math::floorToInt(std::sqrt((Float)m_lobeSize.y));
		m_woResolution = math::floorToInt(std::sqrt((Float)m_lobeSize.x));

		Log(EInfo, "wiRes = %d, woRes = %d", m_wiResolution, m_woResolution);

		BSDF::configure();
	}

	Spectrum evalScale(const BSDFSamplingRecord &bRec) const {
		Vector wiWorld = bRec.its.toWorld(bRec.wi);
		Vector wiMacro = bRec.its.baseFrame.toLocal(wiWorld);
		Vector woWorld = bRec.its.toWorld(bRec.wo);
		Vector woMacro = bRec.its.baseFrame.toLocal(woWorld);

		int rWi = math::floorToInt((wiMacro.y + 1.0) * 0.5 * m_wiResolution);
		rWi = math::clamp(m_wiResolution - rWi - 1, 0, m_wiResolution - 1);
		int cWi = math::floorToInt((wiMacro.x + 1.0) * 0.5 * m_wiResolution);
		cWi = math::clamp(cWi, 0, m_wiResolution - 1);

		int rWo = math::floorToInt((woMacro.y + 1.0) * 0.5 * m_woResolution);
		rWo = math::clamp(m_woResolution - rWo - 1, 0, m_woResolution - 1);
		int cWo = math::floorToInt((woMacro.x + 1.0) * 0.5 * m_woResolution);
		cWo = math::clamp(cWo, 0, m_woResolution - 1);

		Spectrum res(0.f);
		Float weights = 0.f;

		res = m_angularScales->getPixel(Point2i(rWo * m_woResolution + cWo, 
			rWi * m_wiResolution + cWi));

		if (res[0] < 1e-4f) {
// 			Log(EInfo, "(%d, %d), wi = (%.6f, %.6f, %.6f), wo = (%.6f, %.6f, %.6f), (%.6f, %.6f, %.6f)", 
// 				rWi * m_wiResolution + cWi, rWo * m_woResolution + cWo, 
// 				wiMacro.x, wiMacro.y, wiMacro.z,
// 				woMacro.x, woMacro.y, woMacro.z,
// 				res[0], res[1], res[2]);

			res = Spectrum(1.0f);
		}
		//Log(EInfo, "(%.6f, %.6f, %.6f), (%d, %d), %.6f", 
		//	woMacro.x, woMacro.y, woMacro.z, r, c, res.average());
		return res;
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		Spectrum scale = evalScale(bRec);
		if (scale.isZero())
			return Spectrum(0.f);
		Spectrum spec = m_bsdf->eval(bRec, measure);
		return spec * scale;
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		return m_bsdf->pdf(bRec, measure);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		Spectrum spec = m_bsdf->sample(bRec, pdf, sample);
		if (spec.isZero())
			return Spectrum(0.f);
		Spectrum scale = evalScale(bRec);
		return spec * scale;
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sp) const {
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

	std::string toString() const {
		std::ostringstream oss;
		oss << "TabulatedScaledBSDF[" << endl
			<< "  filename = \"" << m_angularScaleFilename << "\"," << endl;
		oss << "]";
		return oss.str();
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

	MTS_DECLARE_CLASS()
public:
	ref<Bitmap> m_angularScales;
	ref<BSDF> m_bsdf;
	std::string m_angularScaleFilename;
	Vector2i m_lobeSize;
	int m_wiResolution;
	int m_woResolution;
};

MTS_IMPLEMENT_CLASS_S(TabulatedScaledBSDF, false, BSDF)
MTS_EXPORT_PLUGIN(TabulatedScaledBSDF, "Tabulated scaled BSDF");
MTS_NAMESPACE_END
