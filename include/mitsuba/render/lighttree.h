/*
	Added by Lifan Wu
	Oct 18, 2015
*/

#pragma once
#if !defined(__MITSUBA_RENDER_LIGHT_TREE_H_)
#define __MITSUBA_RENDER_LIGHT_TREE_H_	

#include <mitsuba/render/lightnode.h>

MTS_NAMESPACE_BEGIN

template <typename LightNodeType> class LightTree {
public:
	LightTree() : root(NULL), 
		nMaxNodes(0), nextFreeNode(1), m_ratio(0), m_random(NULL) {}

	inline void init(Random *random, Float ratio) {
		m_random = random;
		m_ratio = ratio;
	}

	void build(std::vector<VPL> &_vpls);

	LightNodeType *recursiveBuild(size_t nodeNum, int start, int end, 
		VPL **vpls);

	LightNodeType *root;

	size_t nMaxNodes, nextFreeNode;
	std::vector<LightNodeType> m_nodes;

	Random *m_random;

	// control weights between positon and direction, useful in surfaceLight
	Float m_ratio;
};

template <typename LightNodeType> 
void LightTree<LightNodeType>::build(std::vector<VPL> &_vpls) {
	nMaxNodes = _vpls.size() * 4;
	nextFreeNode = 1;
	m_nodes.resize(nMaxNodes);
	std::vector<VPL*> vpls(_vpls.size(), NULL);
	for (size_t i = 0; i < _vpls.size(); i++) {
		vpls[i] = &_vpls[i];
		_vpls[i].pos = _vpls[i].getPos(m_ratio);
	}
	root = recursiveBuild(0, 0, _vpls.size(), &vpls[0]);
}

template <typename LightNodeType>
LightNodeType* LightTree<LightNodeType>::recursiveBuild(size_t nodeNum, int start, int end, 
		VPL **vpls) {
	if (start + 1 == end) {
		m_nodes[nodeNum].initLeaf(vpls[start]);
		return &m_nodes[nodeNum];
	}

	AABB6 bound;
	for (int i = start; i < end; i++) {
		bound.expandBy(vpls[i]->pos);
	}

	int dim = bound.getLargestAxis();
	int mid = (start + end) / 2;
	std::nth_element(&vpls[start], &vpls[mid], &vpls[end], 
		[dim](const VPL *l1, const VPL *l2)->bool {
			return l1->pos[dim] == l2->pos[dim] ? (l1 < l2) : (l1->pos[dim] < l2->pos[dim]);
	});

	LightNodeType *node = &m_nodes[nodeNum];
	int rc = nextFreeNode++;
	node->right = recursiveBuild(rc, mid, end, vpls);
	int lc = nextFreeNode++;
	node->left = recursiveBuild(lc, start, mid, vpls);
	node->update(m_random);
	return node;
}

typedef LightTree<PointLightNode> PointLightTree;
typedef LightTree<DirectionalLightNode> DirectionalLightTree;
typedef LightTree<SurfaceLightNode> SurfaceLightTree;

MTS_NAMESPACE_END


#endif