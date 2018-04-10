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

		BSDF::configure();
	}

	Spectrum evalScale(const BSDFSamplingRecord &bRec) const {
		Vector woWorld = bRec.its.toWorld(bRec.wo);
		Vector woMacro = bRec.its.baseFrame.toLocal(woWorld);
		int r = math::floorToInt((woMacro.y + 1.0) * 0.5 * m_lobeSize.y);
		r = m_lobeSize.y - r - 1;
		int c = math::floorToInt((woMacro.x + 1.0) * 0.5 * m_lobeSize.x);
		Spectrum res = m_angularScales->getPixel(Point2i(c, r));
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
};

MTS_IMPLEMENT_CLASS_S(TabulatedScaledBSDF, false, BSDF)
MTS_EXPORT_PLUGIN(TabulatedScaledBSDF, "Tabulated scaled BSDF");
MTS_NAMESPACE_END
