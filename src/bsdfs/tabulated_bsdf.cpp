#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/vmf.h>
#include "microfacet.h"

MTS_NAMESPACE_BEGIN

// to be update
// default world frame with n=(0,0,1)

class TabulatedBSDF : public BSDF {
public:
	TabulatedBSDF(const Properties &props) : BSDF(props) {
		m_samplesPerEval = props.getInteger("samplersPerEval", 4);
		m_ndfFilename = props.getString("ndfFilename", "");
	}

	TabulatedBSDF(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		m_samplesPerEval = stream->readInt();
		m_ndfFilename = stream->readString();
		m_bsdf = static_cast<BSDF *>(manager->getInstance(stream));

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeInt(m_samplesPerEval);
		stream->writeString(m_ndfFilename);
		manager->serialize(stream, m_bsdf.get());
	}

	void configure() {
		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide);

		m_usesRayDifferentials = false;

		m_samplers.resize(233);
		m_samplers[0] = static_cast<Sampler *>(PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));
		m_samplers[0]->configure();
		for (int i = 1; i < 233; i++) {
			m_samplers[i] = m_samplers[0]->clone();
		}

		// load ndf
		m_cdf.clear();
		m_normals.clear();
		int cnt = 0;
		ref<Bitmap> ndf = new Bitmap(fs::path(m_ndfFilename));
		Vector2i size = ndf->getSize();
		for (int r = 0; r < size.y; r++) {
			Float y = (r + 0.5) / (Float)size.y * 2.0 - 1.0;
			for (int c = 0; c < size.x; c++) {
				Float x = (c + 0.5) / (Float)size.x * 2.0 - 1.0;
				if (x * x + y * y >= 1.0)
					continue;
				Float z = std::sqrt(1.0 - x * x - y * y);

				Float value = ndf->getPixel(Point2i(c, size.y - r - 1)).average();
				if (value < 1e-8)
					continue;

				value /= z;
				m_normals.push_back(Normal(x, y, z));
				if (cnt == 0) {
					m_cdf.push_back(value);
				}
				else {
					m_cdf.push_back(value + m_cdf[cnt - 1]);
				}
				cnt++;
			}
		}

		Float normFactor = m_cdf[cnt - 1];
		for (int i = 0; i < cnt; i++)
			m_cdf[i] /= normFactor;

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

		//std::ostringstream oss;

		for (int i = 0; i < m_samplesPerEval; i++) {
			int idx = (int)(std::lower_bound(m_cdf.begin(), m_cdf.end(), sampler->next1D()) - m_cdf.begin());
			Assert(idx >= 0 && idx < m_cdf.size());

			Vector norm = m_normals[idx];
			Frame nFrame(norm);
			BSDFSamplingRecord bsdfRec(bRec.its, nFrame.toLocal(wiWorld), nFrame.toLocal(woWorld));
			
			res += m_bsdf->eval(bsdfRec);
		}

		return res / (Float)m_samplesPerEval;
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

	std::string toString() const {
		std::ostringstream oss;
		oss << "TabulatedBSDF[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  samplesPerEval = " << m_samplesPerEval << "," << endl
			<< "  baseBSDF = " << m_bsdf->toString() << "," << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	int m_samplesPerEval;
	std::string m_ndfFilename;
	ref<BSDF> m_bsdf;
	ref_vector<Sampler> m_samplers;

	std::vector<Float> m_cdf;
	std::vector<Normal> m_normals;
};

MTS_IMPLEMENT_CLASS_S(TabulatedBSDF, false, BSDF)
MTS_EXPORT_PLUGIN(TabulatedBSDF, "Base BSDF convolve tabulated NDF");
MTS_NAMESPACE_END