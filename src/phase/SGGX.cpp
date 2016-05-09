/*
	Added by Lifan Wu
	Nov 8, 2015
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/core/plugin.h>
#include "microflake_fiber.h"

MTS_NAMESPACE_BEGIN

class SGGXPhaseFunction : public PhaseFunction {
public:
	enum ESGGXPhaseFunctionType {
		ESpecular = 0x01,
		EDiffuse = 0x02,
		EMixed = 0x03,
		EPhong = 0x04
	};

	SGGXPhaseFunction(const Properties &props) : PhaseFunction(props) {
		std::string typeStr = props.getString("sampleType");
		if (typeStr == "specular")
			m_sampleType = ESpecular;
		else if (typeStr == "diffuse")
			m_sampleType = EDiffuse;
		else if (typeStr == "mixed") {
			m_sampleType = EMixed;
			m_mixedWeight = props.getFloat("mixedWeight", 1.f);
			m_mixedWeight = math::clamp(m_mixedWeight, 0.f, 1.f);
		}
		else
			Log(EError, "Unknown SGGX phase function type. Support specular and diffuse.");

		m_stddev = props.getFloat("stddev", -1.f);
		m_fiberDistr = GaussianFiberDistribution(m_stddev);
	}

	SGGXPhaseFunction(Stream *stream, InstanceManager *manager)
		: PhaseFunction(stream, manager) {
		m_sampleType = (ESGGXPhaseFunctionType)stream->readInt();
		m_stddev = stream->readFloat();
		m_fiberDistr = GaussianFiberDistribution(m_stddev);

		m_mixedWeight = stream->readFloat();

		configure();
	}

	virtual ~SGGXPhaseFunction() { }

	void configure() {
		PhaseFunction::configure();
		m_type = EAnisotropic | ENonSymmetric;

		Properties props("independent");
		m_sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
		m_sampler->configure();

		if (m_stddev != -1.f) {
			Float sigma1 = m_fiberDistr.sigmaT(0.f) * 2.f;
			Float sigma2 = sigma1;
			Float sigma3 = m_fiberDistr.sigmaT(1.f) * 2.f;
			D = Matrix3x3(Vector(sigma1 * sigma1, 0, 0), Vector(0, sigma2 * sigma2, 0), Vector(0, 0, sigma3 * sigma3));
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);
		stream->writeInt((int)m_sampleType);
		stream->writeFloat(m_stddev);
		stream->writeFloat(m_mixedWeight);
	}

	Float evalSingleLobe(const PhaseFunctionSamplingRecord &pRec, int lobeIdx) const {
		Vector wi = pRec.wi;
		Vector wo = pRec.wo;

#ifdef USE_STOC_EVAL
		Spectrum s1 = pRec.mRec.s1;
		Spectrum s2 = pRec.mRec.s2;
#else
		Spectrum s1 = pRec.mRec.s1[lobeIdx];
		Spectrum s2 = pRec.mRec.s2[lobeIdx];
#endif

		Float Sxx = s1[0], Syy = s1[1], Szz = s1[2];
		Float Sxy = s2[0], Sxz = s2[1], Syz = s2[2];

		Float res = 1.f;

		Float sqrSum = Sxx * Sxx + Syy * Syy + Szz * Szz + Sxy * Sxy + Sxz * Sxz + Syz * Syz;
		//if (!(Sxx == 0 && Syy == 0 && Szz == 0 && Sxy == 0 && Sxz == 0 && Syz == 0))
		if (fabsf(sqrSum) < 1e-6f)
			return 0;

		if (m_sampleType == ESpecular) {
			Vector H = wi + wo;
			Float length = H.length();

			if (length == 0)
				return 0.f;

			H /= length;
			res *= 0.25f * ndf(H, Sxx, Syy, Szz, Sxy, Sxz, Syz) / sigma(wi, Sxx, Syy, Szz, Sxy, Sxz, Syz);
		}
		else if (m_sampleType == EDiffuse) {
			Vector wm = sampleVNormal(wi, m_sampler, Sxx, Syy, Szz, Sxy, Sxz, Syz);
			res *= 1.f * INV_PI * std::max(0.f, dot(wo, wm));
		}
		else if (m_sampleType == EMixed) {
			Vector H = wi + wo;
			Float length = H.length();

			if (length == 0)
				return 0.f;

			H /= length;
			Float p1 = 0.25f * ndf(H, Sxx, Syy, Szz, Sxy, Sxz, Syz) / sigma(wi, Sxx, Syy, Szz, Sxy, Sxz, Syz);

			Vector wm = sampleVNormal(wi, m_sampler, Sxx, Syy, Szz, Sxy, Sxz, Syz);
			Float p2 = 1.f * INV_PI * std::max(0.f, dot(wo, wm));

			res *= p1 * m_mixedWeight + p2 * (1.f - m_mixedWeight);
		}
		else
			res = 0;

		return res;
	}

	Float eval(const PhaseFunctionSamplingRecord &pRec) const {
#ifdef USE_STOC_EVAL
		return evalSingleLobe(pRec, 0);
#else
		Float res = 0.f;
		int numLobes = pRec.mRec.numLobes;

		for (int i = 0; i < numLobes; i++) {
			if (pRec.mRec.pdfLobe[i] > 0.f)
				res += pRec.mRec.lobeScales[i] * pRec.mRec.pdfLobe[i] * evalSingleLobe(pRec, i);
		}
		return res;
#endif
	}

	Float eval(const PhaseFunctionSamplingRecord &pRec, Float *weightedF) const {
#ifdef USE_STOC_EVAL
		return evalSingleLobe(pRec, 0);
#else
		Float res = 0.f;
		int numLobes = pRec.mRec.numLobes;

		for (int i = 0; i < numLobes; i++) {
			if (pRec.mRec.pdfLobe[i] > 0.f) {
				weightedF[i] = pRec.mRec.pdfLobe[i] * evalSingleLobe(pRec, i);
				res += pRec.mRec.lobeScales[i] * weightedF[i];
			}
			else {
				weightedF[i] = 0.f;
			}
		}
		return res;
#endif
	}

	inline Float sampleSingleLobe(PhaseFunctionSamplingRecord &pRec, Sampler *sampler, int lobeIdx) const {
		Vector wi = pRec.wi;

#ifdef USE_STOC_EVAL
		Spectrum s1 = pRec.mRec.s1;
		Spectrum s2 = pRec.mRec.s2;
#else
		Spectrum s1 = pRec.mRec.s1[lobeIdx];
		Spectrum s2 = pRec.mRec.s2[lobeIdx];
#endif

		Float Sxx = s1[0], Syy = s1[1], Szz = s1[2];
		Float Sxy = s2[0], Sxz = s2[1], Syz = s2[2];

		Float sqrSum = Sxx * Sxx + Syy * Syy + Szz * Szz + Sxy * Sxy + Sxz * Sxz + Syz * Syz;
		//if (!(Sxx == 0 && Syy == 0 && Szz == 0 && Sxy == 0 && Sxz == 0 && Syz == 0))
		if (fabsf(sqrSum) < 1e-6f)
			return 0;

		Vector wm = sampleVNormal(wi, sampler, Sxx, Syy, Szz, Sxy, Sxz, Syz);

		if (m_sampleType == ESpecular) {
			Vector wo = -wi + 2.f * dot(wm, wi) * wm;
			pRec.wo = normalize(wo);
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
			wo = normalize(wo);
			pRec.wo = wo;
			return 1.f;
		}
		else if (m_sampleType == EMixed) {
			Float flag = sampler->next1D();
			if (flag < m_mixedWeight) {
				// specular
				Vector wo = -wi + 2.f * dot(wm, wi) * wm;
				pRec.wo = normalize(wo);
				return 1.f;
			}
			else {
				// diffuse
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
				wo = normalize(wo);
				pRec.wo = wo;
				return 1.f;
			}
		}
		else
			return 0.f;
	}

	inline Float sample(PhaseFunctionSamplingRecord &pRec, Sampler *sampler) const {
#ifdef USE_STOC_EVAL
		return sampleSingleLobe(pRec, sampler, 0);
#else
		Float cdfs[MAX_SGGX_LOBES];
		int numLobes = pRec.mRec.numLobes;

		for (int i = 0; i < numLobes; i++) {
			if (i == 0)
				cdfs[i] = pRec.mRec.pdfLobe[i] * pRec.mRec.lobeScales[i];
			else
				cdfs[i] = cdfs[i - 1] + pRec.mRec.pdfLobe[i] * pRec.mRec.lobeScales[i];
		}

		if (cdfs[numLobes - 1] < 1e-6f)
			return 0.f;

		Float normFactor = cdfs[numLobes - 1];
		for (int i = 0; i < numLobes; i++) {
			cdfs[i] /= normFactor;
		}

		Float rnd = m_sampler->next1D();
		int lobeIdx = std::lower_bound(cdfs, cdfs + numLobes, rnd) - cdfs;

		Assert(lobeIdx < numLobes);

		return sampleSingleLobe(pRec, sampler, lobeIdx) * normFactor;
#endif
	}

	Float sample(PhaseFunctionSamplingRecord &pRec,
		Float &pdf, Sampler *sampler) const {
		if (sample(pRec, sampler) == 0) {
			pdf = 0; return 0.0f;
		}
		
		// need to be normalized, if using lobe scales
		pdf = eval(pRec);
		Float normFactor = 0.f;
		for (int i = 0; i < pRec.mRec.numLobes; i++) {
			normFactor += pRec.mRec.pdfLobe[i] * pRec.mRec.lobeScales[i];
		}
		pdf /= normFactor;
		
		
		return 1.0f * normFactor;
	}

	Float sample(PhaseFunctionSamplingRecord &pRec,
		Float &pdf, Sampler *sampler, Float *weightedF) const {
		if (sample(pRec, sampler) == 0) {
			pdf = 0; return 0.0f;
		}

		// need to be normalized, if using lobe scales
		pdf = eval(pRec, weightedF);
		Float normFactor = 0.f;
		for (int i = 0; i < pRec.mRec.numLobes; i++) {
			normFactor += pRec.mRec.pdfLobe[i] * pRec.mRec.lobeScales[i];
		}
		pdf /= normFactor;


		return 1.0f * normFactor;
	}

	Float pdf(const PhaseFunctionSamplingRecord &pRec) const {
		Float res = eval(pRec);
		Float normFactor = 0.f;
		for (int i = 0; i < pRec.mRec.numLobes; i++) {
			normFactor += pRec.mRec.pdfLobe[i] * pRec.mRec.lobeScales[i];
		}
		res /= normFactor;
		return res;
	}

	Vector sampleVNormal(const Vector &wi, Sampler *sampler, Float Sxx, Float Syy, Float Szz,
		Float Sxy, Float Sxz, Float Syz) const {
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
		return normalize(wm_kji.x * wk + wm_kji.y * wj + wm_kji.z * wi);
	}

	Float sigma(const Vector &wi, Float Sxx, Float Syy, Float Szz,
		Float Sxy, Float Sxz, Float Syz) const {
		Float sigmaSqr = wi.x * wi.x * Sxx + wi.y * wi.y * Syy + wi.z * wi.z * Szz +
			2.f * (wi.x * wi.y * Sxy + wi.x * wi.z * Sxz + wi.y * wi.z * Syz);
		return (sigmaSqr > 0.f) ? sqrtf(sigmaSqr) : 0.f;
	}

	Float sigmaDir(const Vector &d, Float Sxx, Float Syy, Float Szz,
		Float Sxy, Float Sxz, Float Syz) const {
		return sigma(d, Sxx, Syy, Szz, Sxy, Sxz, Syz);
	}

	Float sigmaDir(const Vector &d, const Spectrum &s1,
		const Spectrum &s2) const {
		return sigma(d, s1[0], s1[1], s1[2], s2[0], s2[1], s2[2]);
	}

	Float sigmaDir(const Vector &d, Spectrum *s1, Spectrum *s2, 
		Float *pdfLobe, int numLobes) const {
		Float res = 0.f;
		for (int i = 0; i < numLobes; i++) {
			if (pdfLobe[i] > 0)
				res += pdfLobe[i] * sigma(d, s1[i][0], s1[i][1], s1[i][2], s2[i][0], s2[i][1], s2[i][2]);
		}
		return res;
	}

	Matrix3x3 getD() const {
		return D;
	}

	Float ndf(const Vector &wm, Float Sxx, Float Syy, Float Szz,
		Float Sxy, Float Sxz, Float Syz) const {
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
		return sigmaDir(0.f);
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
	Sampler *m_sampler;

	GaussianFiberDistribution m_fiberDistr;
	Float m_stddev;
	Matrix3x3 D;

	Float m_mixedWeight;
};


MTS_IMPLEMENT_CLASS_S(SGGXPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(SGGXPhaseFunction, "SGGX phase function");
MTS_NAMESPACE_END