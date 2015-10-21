/*
Added by Lifan Wu
Oct 18, 2015
*/

#include <mitsuba/render/lightcutter.h>

MTS_NAMESPACE_BEGIN

Spectrum Lightcutter::evalLightcut(RayDifferential& ray, RadianceQueryRecord& rRec,
	Random *random, int maxCutSize, Float maxErrorRatio) {
	LightNode *node = m_surfaceLightTree->root;
	const BSDF *bsdf = rRec.its.getBSDF(ray);
	Spectrum contrib = evalNodeIllumination(node, rRec.scene, random, -ray.d, rRec.its, bsdf);
	Spectrum error = evalErrorBound(node, -ray.d, rRec.its, bsdf);

	Spectrum L = contrib;
	std::priority_queue<LightcutHeapNode> q;
	q.push(LightcutHeapNode(node, error, contrib));
	int nLeafNodes = 0;

	while (!q.empty() && q.size() + nLeafNodes < maxCutSize) {
		LightcutHeapNode ln = q.top();
		node = ln.node;
		error = ln.error;
		if (ln.errorValue < L.getLuminance() * maxErrorRatio) {
			break;
		}
		else {
			q.pop();
			if (node->isLeaf()) {
				nLeafNodes++;
				continue;
			}
			else {
				if (node->isLeft) {
					L -= ln.contrib;

					contrib = ln.contrib * node->ratio;
					error = evalErrorBound(node->getLeft(), -ray.d, rRec.its, bsdf);
					L += contrib;
					q.push(LightcutHeapNode(node->getLeft(), error, contrib));

					contrib = evalNodeIllumination(node->getRight(), rRec.scene, random,
						-ray.d, rRec.its, bsdf);
					error = evalErrorBound(node->getRight(), -ray.d, rRec.its, bsdf);
					L += contrib;
					q.push(LightcutHeapNode(node->getRight(), error, contrib));
				}
				else {
					L -= ln.contrib;

					contrib = evalNodeIllumination(node->getLeft(), rRec.scene, random,
						-ray.d, rRec.its, bsdf);
					error = evalErrorBound(node->getLeft(), -ray.d, rRec.its, bsdf);
					L += contrib;
					q.push(LightcutHeapNode(node->getLeft(), error, contrib));

					contrib = ln.contrib * node->ratio;
					error = evalErrorBound(node->getRight(), -ray.d, rRec.its, bsdf);
					L += contrib;
					q.push(LightcutHeapNode(node->getRight(), error, contrib));
				}
			}
		}
	}

	count += 1.f;
	avgCutSize += q.size() + nLeafNodes;

	return L;
}

Spectrum Lightcutter::evalNodeIllumination(LightNode *node, const Scene *scene, Random *random,
	Vector &wi, Intersection &its, const BSDF *bsdf) const {
	Spectrum result;
	if (node->getNodeType() == ESurfaceLightNode) {
		VPL *light = ((SurfaceLightNode*)node)->light;
		DirectSamplingRecord dRec(its);
		Vector wo = normalize(light->its.p - its.p);
		Float d2 = (light->its.p - its.p).lengthSquared();
		// d2 = std::max(m_dLimit, d2);
		Float cosWo = absDot(wo, its.shFrame.n);
		Float cosLight = absDot(-wo, light->its.shFrame.n);
		Float G = cosLight / d2;
		G = std::min(G, m_gLimit);
		BSDFSamplingRecord bRec(its, its.toLocal(wo));
		Spectrum f = bsdf->eval(bRec);
		if (G == 0.f || f.isZero())
			return Spectrum(0.f);
		result = f * G * node->P;
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
	}
	else if (node->getNodeType() == EDirectionalLightNode) {

	}
	else if (node->getNodeType() == EPointLightNode) {

	}
	
	return result;
}

Spectrum Lightcutter::evalErrorBound(LightNode *node, Vector &wi, 
	Intersection &its, const BSDF *bsdf) const {
	Spectrum result(1e8f);
	if (node->getNodeType() == ESurfaceLightNode) {
		Float errorBound = 1.f;

		const AABB &bound = node->getBBox();
		Float dist2 = bound.squaredDistanceTo(its.p);
		
		if (dist2 < 1e-6f)
			return Spectrum(1e8f);
		
		errorBound /= dist2;
		
		// Bound cosine term
		if (m_useCosBound) {
			Float cosBound = 1.f;
			Float cosHalfAngle = node->getCone().getAngleCos();

			if (cosHalfAngle > 0.f) {
				Vector coneAxis = node->getCone().getAxis();
				Vector v = cross(Vector(0.f, 0.f, 1.f), coneAxis);
				if (v.length() > 0) {
					v = normalize(v);
					Transform t = Transform::rotate(v, -acosf(dot(Vector(0.f, 0.f, 1.f), coneAxis)));

					Point corners[8];
					AABB tBound;
					Point tp;
					for (int i = 0; i < 8; i++) {
						corners[i] = node->getBBox().getCorner(i);
						t((Point)(its.p - corners[i]), tp);
						tBound.expandBy(tp);
					}

					Point &m = tBound.min;
					Point &M = tBound.max;

					Float cosTheta = 0.f;
					if (M.z > 0) {
						Float minx2, miny2;
						if (m.x * M.x <= 0)
							minx2 = 0.f;
						else
							minx2 = std::min(m.x * m.x, M.x * M.x);

						if (m.y * M.y <= 0)
							miny2 = 0.f;
						else
							miny2 = std::min(m.y * m.y, M.y * M.y);
						Float maxz2 = M.z * M.z;
						cosTheta = M.z / sqrt(minx2 + miny2 + maxz2);
					}

					cosTheta = math::clamp(cosTheta, 0.f, 1.f);

					if (cosTheta <= cosHalfAngle) {
						Float sinHalfAngle = sqrt(1.f - cosHalfAngle * cosHalfAngle);
						Float sinTheta = sqrt(1.f - cosTheta * cosTheta);
						cosBound = math::clamp(cosTheta * cosHalfAngle + sinTheta * sinHalfAngle, 0.f, 1.f);
					}
				}
			}
			errorBound *= cosBound;
		}
		
		Spectrum brdf = evalMaterialErrorBound(node->getBBox(), wi, its, bsdf, false);

		result = node->P * brdf * errorBound;
	}
	else if (node->getNodeType() == EDirectionalLightNode) {

	}
	else if (node->getNodeType() == EPointLightNode) {

	}

	return result;
}

/* 
	Lambertian material
*/
Spectrum Lightcutter::evalMaterialErrorBound(const AABB &bbox, Vector &wi,
	Intersection &its, const BSDF *bsdf, bool isDirLight) const {
	Spectrum diffuse = bsdf->getDiffuseReflectance(its) * INV_PI;
	
	if (!m_useCosBound)
		return diffuse;

	Vector v = cross(Vector(0.f, 0.f, 1.f), its.geoFrame.n);
	if (v.length() < 1e-6f)
		return diffuse;

	v = normalize(v);
	Transform t = Transform::rotate(v, -acosf(dot(Vector(0.f, 0.f, 1.f), its.geoFrame.n)));
	Point corners[8];
	AABB tBound;
	Point tp;
	for (int i = 0; i < 8; i++) {
		corners[i] = bbox.getCorner(i);
		if (!isDirLight) 
			corners[i] = (Point)(corners[i] - its.p);
		else
			corners[i] *= -1;
		t(corners[i], tp);
		tBound.expandBy(tp);
	}

	Point &m = tBound.min;
	Point &M = tBound.max;

	Float cosTheta = 1.f;
	if (M.z > 0) {
		Float minx2, miny2;
		if (m.x * M.x <= 0)
			minx2 = 0.f;
		else
			minx2 = std::min(m.x * m.x, M.x * M.x);

		if (m.y * M.y <= 0)
			miny2 = 0.f;
		else
			miny2 = std::min(m.y * m.y, M.y * M.y);
		Float maxz2 = M.z * M.z;
		cosTheta = M.z / sqrt(minx2 + miny2 + maxz2);
	}
	cosTheta = math::clamp(cosTheta, 0.f, 1.f);

	return diffuse * cosTheta;
}

MTS_NAMESPACE_END