/*
	Added by Lifan Wu
	Nov 8, 2015
*/

#include <mitsuba/core/frame.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/core/plugin.h>
#include "microflake_fiber.h"

MTS_NAMESPACE_BEGIN

class SGGXPhaseFunction : public PhaseFunction {
public:
	enum ESGGXPhaseFunctionType {
		ESpecular = 0x01,
		EDiffuse  = 0x02
	};

	SGGXPhaseFunction(const Properties &props) : PhaseFunction(props) {
		std::string typeStr = props.getString("sampleType");
		if (typeStr == "specular")
			m_sampleType = ESpecular;
		else if (typeStr == "diffuse")
			m_sampleType = EDiffuse;
		else
			Log(EError, "Unknown SGGX phase function type. Support specular and diffuse.");

		Sxx = props.getFloat("Sxx", 0.f);
		Syy = props.getFloat("Syy", 0.f);
		Szz = props.getFloat("Szz", 0.f);
		Sxy = props.getFloat("Sxy", 0.f);
		Sxz = props.getFloat("Sxz", 0.f);
		Syz = props.getFloat("Syz", 0.f);

		m_fiberDistr = GaussianFiberDistribution(0.3f);
	}

	SGGXPhaseFunction(Stream *stream, InstanceManager *manager)
		: PhaseFunction(stream, manager) {
		m_sampleType = (ESGGXPhaseFunctionType)stream->readInt();
		configure();
	}

	virtual ~SGGXPhaseFunction() { }

	void configure() {
		PhaseFunction::configure();
		m_type = EAnisotropic | ENonSymmetric;

		Properties props("independent");
		m_sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
		m_sampler->configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);
		stream->writeInt((int)m_sampleType);
	}

	void setS(Float _Sxx, Float _Syy, Float _Szz, Float _Sxy, Float _Sxz, Float _Syz) {
		Sxx = _Sxx; Syy = _Syy; Szz = _Szz; Sxy = _Sxy; Sxz = _Sxz; Syz = _Syz;
	}

	Float eval(const PhaseFunctionSamplingRecord &pRec) const {
		Vector wi = pRec.wi;
		Vector wo = pRec.wo;

		if (m_sampleType == ESpecular) {
			Vector H = wi + wo;
			Float length = H.length();

			if (length == 0)
				return 0.f;

			H /= length;
			return 0.25f * ndf(H) / sigma(wi);
		}
		else if (m_sampleType == EDiffuse) {
			Vector wm = sampleVNormal(wi, m_sampler);
			return 1.f * INV_PI * std::max(0.f, dot(wo, wm));
		}
		else
			return 0;
	}

	inline Float sample(PhaseFunctionSamplingRecord &pRec, Sampler *sampler) const {
		Vector wi = pRec.wi;
		Vector wm = sampleVNormal(wi, sampler);

		if (m_sampleType == ESpecular) {
			Vector wo = -wi + 2.f * dot(wm, wi) * wm;
			pRec.wo = wo;
			return 1.f;
		}
		else if (m_sampleType == EDiffuse) {
			Float u1 = sampler->next1D();
			Float u2 = sampler->next1D();

			Frame frame(wm);
			Float r1 = 2.f * u1 - 1.f;
			Float r2 = 2.f * u2 - 1.f;

			Float phi, r;
			if (r1 == 0 && r2 == 0) {
				r = phi = 0;
			}
			else if (r1 * r1 > r2 * r2) {
				r = r1;
				phi = (M_PI / 4.f) * (r2 / r1);
			}
			else {
				r = r2;
				phi = (M_PI / 2.f) - (r1 / r2) * (M_PI / 4.f);
			}
			Float x = r * cosf(phi);
			Float y = r * sinf(phi);
			Float z = sqrtf(1.f - x * x - y * y);
			Vector wo = x * frame.s + y * frame.t + z * wm;
			normalize(wo);
			pRec.wo = wo;
			return 1.f;
		}
		else
			return 0.f;
	}

	Float sample(PhaseFunctionSamplingRecord &pRec,
		Float &pdf, Sampler *sampler) const {
		if (sample(pRec, sampler) == 0) {
			pdf = 0; return 0.0f;
		}
		pdf = eval(pRec);
		return 1.0f;
	}

	Vector sampleVNormal(const Vector &wi, Sampler *sampler) const {
		Float u1 = sampler->next1D();
		Float u2 = sampler->next1D();
		
		Float r = sqrtf(u1);
		Float phi = 2.f * M_PI * u2;
		Float u = r * cosf(phi);
		Float v = r * sinf(phi);
		Float w = sqrtf(1.f - u * u - v * v);

		Vector wk, wj;
		Frame frame(wi);
		wk = frame.s;
		wj = frame.t;

		// project S in this basis
		Float Skk = wk.x * wk.x * Sxx + wk.y * wk.y * Syy + wk.z * wk.z * Szz +
			2.f * (wk.x * wk.y * Sxy + wk.x * wk.z * Sxz + wk.y * wk.z * Syz);
		Float Sjj = wj.x * wj.x * Sxx + wj.y * wj.y * Syy + wj.z * wj.z * Szz +
			2.f * (wj.x * wj.y * Sxy + wj.x * wj.z * Sxz + wj.y * wj.z * Syz);
		Float Sii = wi.x * wi.x * Sxx + wi.y * wi.y * Syy + wi.z * wi.z * Szz +
			2.f * (wi.x * wi.y * Sxy + wi.x * wi.z * Sxz + wi.y * wi.z * Syz);
		Float Skj = wk.x * wj.x * Sxx + wk.y * wj.y * Syy + wk.z * wj.z * Szz +
			(wk.x * wj.y + wk.y * wj.x) * Sxy + (wk.x * wj.z + wk.z * wj.x) * Sxz +
			(wk.y * wj.z + wk.z * wj.y) * Syz;
		Float Ski = wk.x * wi.x * Sxx + wk.y * wi.y * Syy + wk.z * wi.z * Szz +
			(wk.x * wi.y + wk.y * wi.x) * Sxy + (wk.x * wi.z + wk.z * wi.x) * Sxz +
			(wk.y * wi.z + wk.z * wi.y) * Syz;
		Float Sji = wj.x * wi.x * Sxx + wj.y * wi.y * Syy + wj.z * wi.z * Szz +
			(wj.x * wi.y + wj.y * wi.x) * Sxy + (wj.x * wi.z + wj.z * wi.x) * Sxz +
			(wj.y * wi.z + wj.z * wi.y) * Syz;

		Float sqrtDetSkji = sqrtf(fabs(Skk * Sjj * Sii - Skj * Skj * Sii - Ski * Ski * Sjj -
			Sji * Sji * Skk + 2.f * Skj * Ski * Sji));
		Float invSqrtSii = 1.f / sqrt(Sii);
		Float tmp = sqrtf(Sjj * Sii - Sji * Sji);
		Vector Mk(sqrtDetSkji / tmp, 0.f, 0.f);
		Vector Mj(-invSqrtSii * (Ski * Sji - Skj * Sii) / tmp, invSqrtSii * tmp, 0);
		Vector Mi(invSqrtSii * Ski, invSqrtSii * Sji, invSqrtSii * Sii);
		
		Vector wm_kji = normalize(u * Mk + v * Mj + w * Mi);
		return wm_kji.x * wk + wm_kji.y * wj + wm_kji.z * wi;
	}

	Float sigma(const Vector &wi) const {
		Float sigmaSqr = wi.x * wi.x * Sxx + wi.y * wi.y * Syy + wi.z * wi.z * Szz +
			2.f * (wi.x * wi.y * Sxy + wi.x * wi.z * Sxz + wi.y * wi.z * Syz);
		return (sigmaSqr > 0.f) ? sqrtf(sigmaSqr) : 0.f;
	}

	Float ndf(const Vector &wm) const {
		Float detS = Sxx * Syy * Szz - Sxx * Syz * Syz - Syy * Sxz * Sxz - Szz * Sxy * Sxy + 2.f * Sxy * Sxz * Syz;
		Float den = wm.x * wm.x * (Syy * Szz - Syz * Syz) + wm.y * wm.y * (Sxx * Szz - Sxz * Sxz) +
			wm.z * wm.z * (Sxx * Syy - Sxy * Sxy) + 2.f * (wm.x * wm.y * (Sxz * Syz - Szz * Sxy) +
			wm.x * wm.z * (Sxy * Syz - Syy * Sxz) + wm.y * wm.z * (Sxy * Sxz - Sxx * Syz));
		return powf(fabsf(detS), 1.5f) / (M_PI * den * den);
	}

	bool needsDirectionallyVaryingCoefficients() const { return true; }

	Float sigmaDir(Float cosTheta) const {
		// Scaled such that replacing an isotropic phase function with an
		// isotropic microflake distribution does not cause changes
		return 2 * m_fiberDistr.sigmaT(cosTheta);
	}

	Float sigmaDirMax() const {
		return sigmaDir(0);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MicroflakePhaseFunction[" << endl
			<< "   sampleType = " << m_sampleType << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ESGGXPhaseFunctionType m_sampleType;
	Float Sxx, Syy, Szz, Sxy, Sxz, Syz;
	Sampler *m_sampler;

	GaussianFiberDistribution m_fiberDistr;
};


MTS_IMPLEMENT_CLASS_S(SGGXPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(SGGXPhaseFunction, "SGGX phase function");
MTS_NAMESPACE_END