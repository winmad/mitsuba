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

	// Luminance attenuation of a child
	Float ratio;
	bool isLeft;
	
	virtual LightNode* getLeft() const = 0;
	virtual LightNode* getRight() const = 0;
	virtual AABB getBBox() const = 0;
	virtual VPL* getLight() const = 0;
	virtual Cone getCone() const = 0;

	virtual bool isLeaf() const = 0;
	virtual ELightNodeType getNodeType() const = 0;
	virtual void initLeaf(VPL *vl) = 0;
	virtual void update(Random *random) = 0;
};

struct PointLightNode : public LightNode {
	PointLightNode() : LightNode(), light(NULL), left(NULL), right(NULL) {}
	PointLightNode(VPL &vl) : light(&vl), bbox(vl.its.p) {
		P = vl.P;
	}
	virtual ~PointLightNode() {}
	
	virtual LightNode* getLeft() const {
		return left;
	}

	virtual LightNode* getRight() const {
		return right;
	}

	virtual AABB getBBox() const {
		return bbox;
	}

	virtual VPL* getLight() const {
		return light;
	}

	virtual Cone getCone() const {
		return Cone();
	}

	virtual bool isLeaf() const {
		return left == NULL && right == NULL;
	}

	virtual ELightNodeType getNodeType() const {
		return EPointLightNode;
	}

	virtual void initLeaf(VPL *vl) {
		left = right = NULL;
		light = vl;
		P = vl->P;
		bbox.reset();
		bbox.expandBy(vl->its.p);
	}

	virtual void update(Random *random) {
		bbox = left->bbox;
		bbox.expandBy(right->bbox);
		P = left->P + right->P;
		ratio = left->P.getLuminance() / P.getLuminance();
		if (random->nextFloat() < ratio) {
			light = left->light;
			isLeft = true;
		}
		else {
			light = right->light;
			ratio = 1.f - ratio;
			isLeft = false;
		}
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

	virtual LightNode* getLeft() const {
		return left;
	}

	virtual LightNode* getRight() const {
		return right;
	}

	virtual AABB getBBox() const {
		return bbox;
	}

	virtual VPL* getLight() const {
		return light;
	}

	virtual Cone getCone() const {
		return Cone();
	}

	virtual bool isLeaf() const {
		return left == NULL && right == NULL;
	}

	virtual ELightNodeType getNodeType() const {
		return EDirectionalLightNode;
	}

	virtual void initLeaf(VPL *vl) {
		left = right = NULL;
		light = vl;
		P = vl->P;
		bbox.reset();
		bbox.expandBy(vl->its.p);
	}

	virtual void update(Random *random) {
		bbox = left->bbox;
		bbox.expandBy(right->bbox);
		P = left->P + right->P;
		ratio = left->P.getLuminance() / P.getLuminance();
		if (random->nextFloat() < ratio) {
			light = left->light;
			isLeft = true;
		}
		else {
			light = right->light;
			ratio = 1.f - ratio;
			isLeft = false;
		}
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

	virtual LightNode* getLeft() const {
		return left;
	}

	virtual LightNode* getRight() const {
		return right;
	}

	virtual AABB getBBox() const {
		return bbox;
	}

	virtual VPL* getLight() const {
		return light;
	}

	virtual Cone getCone() const {
		return cone;
	}

	virtual bool isLeaf() const {
		return left == NULL && right == NULL;
	}

	virtual ELightNodeType getNodeType() const {
		return ESurfaceLightNode;
	}

	virtual void initLeaf(VPL *vl) {
		left = right = NULL;
		light = vl;
		
		// Assume Lambertian BRDF
		Spectrum f;
		if (vl->type == ESurfaceVPL)
			f = vl->its.getBSDF()->getDiffuseReflectance(vl->its) * INV_PI;
		else
			f = Spectrum(INV_PI);
		
		P = vl->P * f;
		bbox.reset();
		bbox.expandBy(vl->its.p);
		cone.reset();
		cone.expandBy(vl->its.shFrame.n);
	}

	virtual void update(Random *random) {
		bbox = left->bbox;
		bbox.expandBy(right->bbox);
		cone = left->cone;
		cone.expandBy(right->cone);
		P = left->P + right->P;
		ratio = left->P.getLuminance() / P.getLuminance();
		if (random->nextFloat() < ratio) {
			light = left->light;
			isLeft = true;
		}
		else {
			light = right->light;
			ratio = 1.f - ratio;
			isLeft = false;
		}
	}

	AABB bbox;
	Cone cone;
	VPL *light;
	SurfaceLightNode *left, *right;
};

MTS_NAMESPACE_END

#endif