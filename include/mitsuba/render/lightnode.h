/*
	Added by Lifan Wu
	Oct 17, 2015
*/

#pragma once
#if !defined(__MITSUBA_RENDER_LIGHT_NODE_H_)
#define __MITSUBA_RENDER_LIGHT_NODE_H_	

#include <mitsuba/render/vpl.h>
#include <mitsuba/core/cone.h>

MTS_NAMESPACE_BEGIN

enum ELightNodeType {
	EPointLightNode = 0,
	EDirectionalLightNode,
	ESurfaceLightNode
};

struct LightNode {
	LightNode() : P(0.f) {}
	virtual ~LightNode() {}
	
	Spectrum P;
	
	virtual bool isLeaf() const = 0;
	virtual ELightNodeType getNodeType() const = 0;
	virtual void initLeaf() = 0;
	virtual void update(Random *random) = 0;
};

struct PointLightNode : public LightNode {
	PointLightNode() : LightNode(), light(NULL), left(NULL), right(NULL) {}
	PointLightNode(VPL &vl) : light(&vl), bbox(vl.its.p) {
		P = vl.P;
	}
	virtual ~PointLightNode() {}
	
	virtual bool isLeaf() const {
		return left == NULL && right == NULL;
	}

	virtual ELightNodeType getNodeType() const {
		return EPointLightNode;
	}

	virtual void initLeaf() {
		left = right = NULL;
		light = NULL;
		bbox.reset();
	}

	virtual void update(Random *random) {
		bbox = left->bbox;
		bbox.expandBy(right->bbox);
		P = left->P + right->P;
		light = (random->nextFloat() < (left->P.getLuminance() / P.getLuminance()) ? left->light : right->light);
	}

	AABB bbox;
	VPL *light;
	PointLightNode *left, *right;
};

struct DirectionalLightNode : public LightNode {
	DirectionalLightNode() : LightNode(), light(NULL), left(NULL), right(NULL) {}
	DirectionalLightNode(VPL &vl) : light(&vl), bbox(vl.its.p) {
		P = vl.P;
	}
	virtual ~DirectionalLightNode() {}

	virtual bool isLeaf() const {
		return left == NULL && right == NULL;
	}

	virtual ELightNodeType getNodeType() const {
		return EDirectionalLightNode;
	}

	virtual void initLeaf() {
		left = right = NULL;
		light = NULL;
		bbox.reset();
	}

	virtual void update(Random *random) {
		bbox = left->bbox;
		bbox.expandBy(right->bbox);
		P = left->P + right->P;
		light = (random->nextFloat() < (left->P.getLuminance() / P.getLuminance()) ? left->light : right->light);
	}

	AABB bbox;
	VPL *light;
	DirectionalLightNode *left, *right;
};

struct SurfaceLightNode : public LightNode {
	SurfaceLightNode() : LightNode(), light(NULL), left(NULL), right(NULL) {}
	SurfaceLightNode(VPL &vl) : light(&vl), bbox(vl.its.p), cone(vl.its.shFrame.n) {
		P = vl.P;
	}

	virtual ~SurfaceLightNode() {}

	virtual bool isLeaf() const {
		return left == NULL && right == NULL;
	}

	virtual ELightNodeType getNodeType() const {
		return ESurfaceLightNode;
	}

	virtual void initLeaf() {
		left = right = NULL;
		light = NULL;
		bbox.reset();
		cone.reset();
	}

	virtual void update(Random *random) {
		bbox = left->bbox;
		bbox.expandBy(right->bbox);
		cone = left->cone;
		cone.expandBy(right->cone);
		P = left->P + right->P;
		light = (random->nextFloat() < (left->P.getLuminance() / P.getLuminance()) ? left->light : right->light);
	}

	AABB bbox;
	Cone cone;
	VPL *light;
	SurfaceLightNode *left, *right;
};

MTS_NAMESPACE_END

#endif