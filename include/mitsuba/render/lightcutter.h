/*
Added by Lifan Wu
Oct 18, 2015
*/

#pragma once
#if !defined(__MITSUBA_RENDER_LIGHT_CUTTER_H_)
#define __MITSUBA_RENDER_LIGHT_CUTTER_H_	

#include <mitsuba/render/lighttree.h>

MTS_NAMESPACE_BEGIN

struct LightcutHeapNode {
	LightcutHeapNode() {}
	LightcutHeapNode(LightNode *_node, Spectrum _err, Spectrum _contrib) : node(_node),
		error(_err), contrib(_contrib) {
		errorValue = _err.getLuminance();
	}

	bool operator<(const LightcutHeapNode &d) const {
		if (errorValue == d.errorValue)
			return node > d.node;
		else
			return errorValue < d.errorValue;
	}

	LightNode *node;
	Spectrum error;
	Float errorValue;
	Spectrum contrib;
};

class MTS_EXPORT_RENDER Lightcutter {
public:
	Lightcutter() {}
	Lightcutter(LightTree<PointLightNode> *plt, 
		LightTree<DirectionalLightNode> *dlt, 
		LightTree<SurfaceLightNode> *slt,
		Float gLimit, Float rrThreshold)
		: m_pointLightTree(plt), m_directionalLightTree(dlt),
		m_surfaceLightTree(slt), m_gLimit(gLimit), m_rrThreshold(rrThreshold) {}

	void init(LightTree<PointLightNode> *plt,
		LightTree<DirectionalLightNode> *dlt,
		LightTree<SurfaceLightNode> *slt,
		Float gLimit, Float rrThreshold) {
		m_pointLightTree = plt;
		m_directionalLightTree = dlt;
		m_surfaceLightTree = slt;
		m_gLimit = gLimit;
		m_rrThreshold = rrThreshold;
	}

	Spectrum evalLightcut(const RayDifferential& ray, RadianceQueryRecord& rRec,
		Random *random, int maxCutSize, Float maxErrorRatio) const;
	
	Spectrum evalNodeIllumination(const PointLightNode *node, const Scene *scene, Random *random,
		Vector &wi, Intersection &its, const BSDF *bsdf) const;
	Spectrum evalNodeIllumination(const DirectionalLightNode *node, const Scene *scene, Random *random,
		Vector &wi, Intersection &its, const BSDF *bsdf) const;
	Spectrum evalNodeIllumination(const SurfaceLightNode *node, const Scene *scene, Random *random,
		Vector &wi, Intersection &its, const BSDF *bsdf) const;

	Spectrum evalErrorBound(const PointLightNode *node, Vector &wi,
		Intersection &its, const BSDF *bsdf) const;
	Spectrum evalErrorBound(const DirectionalLightNode *node, Vector &wi,
		Intersection &its, const BSDF *bsdf) const;
	Spectrum evalErrorBound(const SurfaceLightNode *node, Vector &wi,
		Intersection &its, const BSDF *bsdf) const;

	Spectrum evalMaterialErrorBound(const AABB &bbox, Vector &wi,
		Intersection &its, const BSDF *bsdf, bool isDirLight) const;

	LightTree<PointLightNode> *m_pointLightTree;
	LightTree<DirectionalLightNode> *m_directionalLightTree;
	LightTree<SurfaceLightNode> *m_surfaceLightTree;
	Float m_gLimit;
	Float m_rrThreshold;
        
};

MTS_NAMESPACE_END

#endif