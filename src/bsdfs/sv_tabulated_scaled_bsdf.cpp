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

class SVTabulatedScaledBSDF : public BSDF {
public:
	SVTabulatedScaledBSDF(const Properties &props) : BSDF(props) {
		m_xResolution = props.getInteger("xResolution", 1);
		m_yResoluiton = props.getInteger("yResolution", 1);
		Float uvscale = props.getFloat("uvscale", 1.0f);
		m_uvScale = Vector2(
			props.getFloat("uscale", uvscale),
			props.getFloat("vscale", uvscale)
			);
		m_angularScalePrefix = props.getString("angularScalePrefix", "");
	}

	SVTabulatedScaledBSDF(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		m_xResolution = stream->readInt();
		m_yResoluiton = stream->readInt();
		m_uvScale = Vector2(stream->readFloat(), stream->readFloat());
		m_angularScalePrefix = stream->readString();

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);
		stream->writeInt(m_xResolution);
		stream->writeInt(m_yResoluiton);
		stream->writeFloat(m_uvScale.x);
		stream->writeFloat(m_uvScale.y);
		stream->writeString(m_angularScalePrefix);
	}

	void configure() {
		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide | ESpatiallyVarying);

		m_usesRayDifferentials = false;

		// load angular scales
		m_angularScales.resize(m_xResolution * m_yResoluiton);
		for (int y = 0; y < m_yResoluiton; y++) {
			for (int x = 0; x < m_xResolution; x++) {
				std::ostringstream oss;
				oss << m_angularScalePrefix << "_" << x << "_" << y << ".exr";
				m_angularScales[y * m_xResolution + x] = new Bitmap(fs::path(oss.str()));
			}
		}
		Log(EInfo, "sv scaling matrix loaded.");

		m_lobeSize = m_angularScales[0]->getSize();
		m_wiResolution = math::floorToInt(std::sqrt((Float)m_lobeSize.y));
		m_woResolution = math::floorToInt(std::sqrt((Float)m_lobeSize.x));
		Log(EInfo, "wiRes = %d, woRes = %d", m_wiResolution, m_woResolution);

		BSDF::configure();
	}

	Spectrum evalScale(const BSDFSamplingRecord &bRec) const {
		Vector2 uv;
		uv[0] = bRec.its.uv[0] * m_uvScale[0];
		uv[1] = bRec.its.uv[1] * m_uvScale[1];
		uv[0] = uv[0] - math::floorToInt(uv[0]);
		uv[1] = uv[1] - math::floorToInt(uv[1]);
		int scaleMatIdx = math::floorToInt(uv[1] * m_yResoluiton) * m_xResolution 
			+ math::floorToInt(uv[0] * m_xResolution);

		Vector wiWorld = bRec.its.toWorld(bRec.wi);
		Vector wiMacro = bRec.its.baseFrame.toLocal(wiWorld);
		Vector woWorld = bRec.its.toWorld(bRec.wo);
		Vector woMacro = bRec.its.baseFrame.toLocal(woWorld);

		if (wiMacro.z <= 0 || woMacro.z <= 0)
			return Spectrum(0.0);

		Point2 wiTex = warp::uniformHemisphereToSquareConcentric(wiMacro);
		Point2 woTex = warp::uniformHemisphereToSquareConcentric(woMacro);

		int wiNumCells = m_wiResolution - 1;
		int woNumCells = m_woResolution - 1;

		int c1 = math::clamp(math::floorToInt(wiTex.x * wiNumCells), 0, wiNumCells - 1);
		int r1 = math::clamp(math::floorToInt(wiTex.y * wiNumCells), 0, wiNumCells - 1);
		int c2 = math::clamp(math::floorToInt(woTex.x * woNumCells), 0, woNumCells - 1);
		int r2 = math::clamp(math::floorToInt(woTex.y * woNumCells), 0, woNumCells - 1);

		Spectrum res(0.f);
		for (int dr1 = 0; dr1 < 2; dr1++) {
			Float v1 = wiTex.y * wiNumCells - r1;
			Float wv1 = std::abs(1.0 - dr1 - v1);

			for (int dc1 = 0; dc1 < 2; dc1++) {
				Float u1 = wiTex.x * wiNumCells - c1;
				Float wu1 = std::abs(1.0 - dc1 - u1);

				for (int dr2 = 0; dr2 < 2; dr2++) {
					Float v2 = woTex.y * woNumCells - r2;
					Float wv2 = std::abs(1.0 - dr2 - v2);

					for (int dc2 = 0; dc2 < 2; dc2++) {
						Float u2 = woTex.x * woNumCells - c2;
						Float wu2 = std::abs(1.0 - dc2 - u2);

						int wiIdx = (r1 + dr1) * m_wiResolution + (c1 + dc1);
						int woIdx = (r2 + dr2) * m_woResolution + (c2 + dc2);

						Spectrum tmpValue = m_angularScales[scaleMatIdx]->getPixel(Point2i(woIdx, wiIdx));
						res += tmpValue * wv1 * wu1 * wv2 * wu2;
					}
				}
			}
		}

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
		oss << "SVTabulatedScaledBSDF[" << endl
			<< "  filenamePrefix = \"" << m_angularScalePrefix << "\"," << endl;
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
	int m_xResolution;
	int m_yResoluiton;
	Vector2 m_uvScale;
	
	std::string m_angularScalePrefix;
	ref_vector<Bitmap> m_angularScales;
	
	ref<BSDF> m_bsdf;
	Vector2i m_lobeSize;
	int m_wiResolution;
	int m_woResolution;
};

MTS_IMPLEMENT_CLASS_S(SVTabulatedScaledBSDF, false, BSDF)
MTS_EXPORT_PLUGIN(SVTabulatedScaledBSDF, "Spatially-varying tabulated scaled BSDF");
MTS_NAMESPACE_END
