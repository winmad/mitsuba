#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN

std::string PhaseFunctionSamplingRecord::toString() const {
	std::ostringstream oss;
	oss << "PhaseFunctionSamplingRecord[" << endl
		<< "  mRec = " << indent(mRec.toString()) << "," << endl
		<< "  wi = " << wi.toString() << "," << endl
		<< "  wo = " << wo.toString() << "," << endl
		<< "  mode = " << mode << endl
		<< "]";
	return oss.str();
}

PhaseFunctionSamplingRecord::PhaseFunctionSamplingRecord(const MediumSamplingRecord &mRec,
	const Vector &wi, bool _useSGGX, ETransportMode mode)
	: mRec(mRec), wi(wi), mode(mode), useSGGX(_useSGGX) {
	if (useSGGX) {
		const VolumeDataSource *VolS1 = mRec.medium->getS1();
		const VolumeDataSource *VolS2 = mRec.medium->getS2();
		if (VolS1 == NULL) {
			Matrix3x3 D = mRec.medium->getPhaseFunction()->getD();
			Vector w3 = mRec.orientation;
			Frame frame(w3);
			Matrix3x3 basis(frame.s, frame.t, w3);
			Matrix3x3 basisT;
			basis.transpose(basisT);
			Matrix3x3 S = basis * D * basisT;
			Sxx = S.m[0][0]; Syy = S.m[1][1]; Szz = S.m[2][2];
			Sxy = S.m[0][1]; Sxz = S.m[0][2]; Syz = S.m[1][2];
		}
		else {
			Spectrum S1 = VolS1->lookupSpectrum(mRec.p);
			Spectrum S2 = VolS2->lookupSpectrum(mRec.p);
			// fix SGGX volumeToWorld rotation!!!
			Sxx = S1[0]; Syy = S1[1]; Szz = S1[2];
			Sxy = S2[0]; Sxz = S2[1]; Syz = S2[2];
		}
	}
}

PhaseFunctionSamplingRecord::PhaseFunctionSamplingRecord(const MediumSamplingRecord &mRec,
	const Vector &wi, const Vector &wo, bool _useSGGX, ETransportMode mode)
	: mRec(mRec), wi(wi), wo(wo), mode(mode), useSGGX(_useSGGX) {
	if (useSGGX) {
		const VolumeDataSource *VolS1 = mRec.medium->getS1();
		const VolumeDataSource *VolS2 = mRec.medium->getS2();
		if (VolS1 == NULL) {
			Matrix3x3 D = mRec.medium->getPhaseFunction()->getD();
			Vector w3 = mRec.orientation;
			Frame frame(w3);
			Matrix3x3 basis(frame.s, frame.t, w3);
			Matrix3x3 basisT;
			basis.transpose(basisT);
			Matrix3x3 S = basis * D * basisT;
			Sxx = S.m[0][0]; Syy = S.m[1][1]; Szz = S.m[2][2];
			Sxy = S.m[0][1]; Sxz = S.m[0][2]; Syz = S.m[1][2];
		}
		else {
			Spectrum S1 = VolS1->lookupSpectrum(mRec.p);
			Spectrum S2 = VolS2->lookupSpectrum(mRec.p);
			// fix SGGX volumeToWorld rotation!!!
			Sxx = S1[0]; Syy = S1[1]; Szz = S1[2];
			Sxy = S2[0]; Sxz = S2[1]; Syz = S2[2];
		}
	}
}

void PhaseFunction::configure() {
	m_type = 0;
}

Float PhaseFunction::pdf(const PhaseFunctionSamplingRecord &pRec) const {
	return eval(pRec);
}

bool PhaseFunction::needsDirectionallyVaryingCoefficients() const {
	return false;
}

Float PhaseFunction::sigmaDir(Float cosTheta) const {
	Log(EError, "%s::sigmaDir(Float) is not implemented (this is not "
		"an anisotropic medium!)", getClass()->getName().c_str());
	return 0.0f;
}

Matrix3x3 PhaseFunction::getD() const {
	Log(EError, "%s::getD() is not implemented (Only for SGGX)", 
		getClass()->getName().c_str());
	return Matrix3x3();
}

Float PhaseFunction::sigmaDir(const Vector &d, Float Sxx, Float Syy, Float Szz,
	Float Sxy, Float Sxz, Float Syz) const {
	Log(EError, "%s::sigmaDir(Vector) is not implemented (this is not "
		"an anisotropic medium!)", getClass()->getName().c_str());
	return 0.0f;
}

Float PhaseFunction::sigmaDirMax() const {
	Log(EError, "%s::sigmaDirMax() is not implemented (this is not "
		"an anisotropic medium!)", getClass()->getName().c_str());
	return 0.0f;
}

Float PhaseFunction::getMeanCosine() const {
	Log(EError, "%s::getMeanCosine() is not implemented!",
		getClass()->getName().c_str());
	return 0.0f;
}

MTS_IMPLEMENT_CLASS(PhaseFunction, true, ConfigurableObject)
MTS_NAMESPACE_END
