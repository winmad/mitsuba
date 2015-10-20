/*
Added by Lifan Wu
Oct 18, 2015
*/

#include <mitsuba/render/lightcutter.h>

MTS_NAMESPACE_BEGIN

Spectrum Lightcutter::evalLightcut(const RayDifferential& ray, RadianceQueryRecord& rRec,
	Random *random, int maxCutSize, Float maxErrorRatio) const {
	LightNode *node = m_surfaceLightTree->root;
	const BSDF *bsdf = rRec.its.getBSDF(ray);
	Spectrum contrib = evalNodeIllumination((SurfaceLightNode*)node, rRec.scene, random, -ray.d, rRec.its, bsdf);
	Spectrum error = evalErrorBound((SurfaceLightNode*)node, -ray.d, rRec.its, bsdf);

	Spectrum L = contrib;
	std::priority_queue<LightcutHeapNode> q;
	q.push(LightcutHeapNode(node, error, contrib));
	int nLeafNodes = 0;

	while (!q.empty() && q.size() + nLeafNodes < maxCutSize) {
		LightcutHeapNode ln = q.top();
		node = ln.node;
		error = ln.error;
		if (error[0] < L[0] * maxErrorRatio && error[1] < L[1] * maxErrorRatio &&
			error[2] < L[2] * maxErrorRatio && error.getLuminance() < L.getLuminance() * maxErrorRatio) {
			break;
		}
		else {
			q.pop();
			if (node->isLeaf()) {
				nLeafNodes++;
				continue;
			}
			else {
				L -= ln.contrib;
				
				SurfaceLightNode *left = ((SurfaceLightNode*)node)->left;
				contrib = evalNodeIllumination(left, rRec.scene, random,
					-ray.d, rRec.its, bsdf);
				error = evalErrorBound(left, -ray.d, rRec.its, bsdf);
				L += contrib;
				q.push(LightcutHeapNode(left, error, contrib));

				SurfaceLightNode *right = ((SurfaceLightNode*)node)->right;
				contrib = evalNodeIllumination(right, rRec.scene, random,
					-ray.d, rRec.its, bsdf);
				error = evalErrorBound(right, -ray.d, rRec.its, bsdf);
				L += contrib;
				q.push(LightcutHeapNode(right, error, contrib));
			}
		}
	}

	return L;
}

Spectrum Lightcutter::evalNodeIllumination(const PointLightNode *node, const Scene *scene, Random *random,
	Vector &wi, Intersection &its, const BSDF *bsdf) const {
	return Spectrum(0.f);
}

Spectrum Lightcutter::evalNodeIllumination(const DirectionalLightNode *node, const Scene *scene, Random *random,
	Vector &wi, Intersection &its, const BSDF *bsdf) const {
	return Spectrum(0.f);
}

Spectrum Lightcutter::evalNodeIllumination(const SurfaceLightNode *node, const Scene *scene, Random *random,
	Vector &wi, Intersection &its, const BSDF *bsdf) const {
	VPL *light = node->light;
	DirectSamplingRecord dRec(its);
	Vector wo = normalize(light->its.p - its.p);
	Float d2 = (light->its.p - its.p).lengthSquared();
	Float cosWo = absDot(wo, its.shFrame.n);
	Float cosLight = absDot(-wo, light->its.shFrame.n);
	Float G = cosLight / d2;
	G = std::min(G, m_gLimit);
	BSDFSamplingRecord bRec(its, its.toLocal(wo));
	Spectrum f = bsdf->eval(bRec);
	if (G == 0.f || f.isZero()) 
		return Spectrum(0.f);
	Spectrum result = f * G * node->P;
	Ray connectRay(its.p, wo, Epsilon, sqrt(d2) * (1.f - ShadowEpsilon), its.time);

	if (result.getLuminance() < m_rrThreshold) {
		Float continueProb = 0.1f;
		if (random->nextFloat() > continueProb)
			return Spectrum(0.f);
		result /= continueProb;
	}
	
	if (!scene->rayIntersect(connectRay))
		return result;
	else
		return Spectrum(0.f);
	
	return result;
}

Spectrum Lightcutter::evalErrorBound(const SurfaceLightNode *node, Vector &wi, 
	Intersection &its, const BSDF *bsdf) const {
	Float errorBound = 1.f;

	const AABB &bound = node->bbox;
	Float dist2 = bound.squaredDistanceTo(its.p);
	if (dist2 < 1e-5f)
		return Spectrum(1e8f);
	errorBound /= dist2;

	return node->P * errorBound;
}

MTS_NAMESPACE_END