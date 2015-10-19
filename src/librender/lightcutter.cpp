/*
Added by Lifan Wu
Oct 18, 2015
*/

#include <mitsuba/render/lightcutter.h>

MTS_NAMESPACE_BEGIN

Spectrum Lightcutter::evalLightcut(const RayDifferential& ray, RadianceQueryRecord& rRec,
	int maxCutSize, Float maxErrorRatio) {
	
}

Spectrum Lightcutter::evalNodeIllumination(const SurfaceLightNode *node, Vector &wi,
	Intersection &its, const BSDF *bsdf) {
	return Spectrum(0.f);
}

Spectrum Lightcutter::evalNodeIllumination(const DirectionalLightNode *node, Vector &wi,
	Intersection &its, const BSDF *bsdf) {
	return Spectrum(0.f);
}

Spectrum Lightcutter::evalNodeIllumination(const SurfaceLightNode *node, Vector &wi,
	Intersection &its, const BSDF *bsdf) {

}

MTS_NAMESPACE_END