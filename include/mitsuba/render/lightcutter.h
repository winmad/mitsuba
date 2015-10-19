/*
Added by Lifan Wu
Oct 18, 2015
*/

#pragma once
#if !defined(__MITSUBA_RENDER_LIGHT_CUTTER_H_)
#define __MITSUBA_RENDER_LIGHT_CUTTER_H_	

#include <mitsuba/render/lighttree.h>

MTS_NAMESPACE_BEGIN

class Lightcutter {
public:
	Lightcutter(const LightTree<PointLightNode> *plt, 
		const LightTree<DirectionalLightNode> *dlt, 
		const LightTree<SurfaceLightNode> *slt,
		Float gLimit = 0.01f, Float rrThreshold = 1e-4f)
		: m_pointLightTree(plt), m_directionalLightTree(dlt),
		m_surfaceLightTree(slt), m_gLimit(gLimit), m_rrThreshold(rrThreshold) {}

	Spectrum evalLightcut(const RayDifferential& ray, RadianceQueryRecord& rRec,
		int maxCutSize, Float maxErrorRatio);
	
	Spectrum evalNodeIllumination(const PointLightNode *node, Vector &wi, 
		Intersection &its, const BSDF *bsdf);
	Spectrum evalNodeIllumination(const DirectionalLightNode *node, Vector &wi,
		Intersection &its, const BSDF *bsdf);
	Spectrum evalNodeIllumination(const SurfaceLightNode *node, Vector &wi,
		Intersection &its, const BSDF *bsdf);

	Spectrum evalErrorBound(const PointLightNode *node, Vector &wi,
		Intersection &its, const BSDF *bsdf);
	Spectrum evalErrorBound(const DirectionalLightNode *node, Vector &wi,
		Intersection &its, const BSDF *bsdf);
	Spectrum evalErrorBound(const SurfaceLightNode *node, Vector &wi,
		Intersection &its, const BSDF *bsdf);

	Spectrum evalMaterialErrorBound(const AABB &bbox, Vector &wi,
		Intersection &its, const BSDF *bsdf, bool isDirLight);

	const LightTree<PointLightNode> *m_pointLightTree;
	const LightTree<DirectionalLightNode> *m_directionalLightTree;
	const LightTree<SurfaceLightNode> *m_surfaceLightTree;
	Float m_gLimit;
	Float m_rrThreshold;
};

MTS_NAMESPACE_END

#endif