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
	LightTree() : root(NULL), nNodes(0), nextFreeNode(1) {}

	void build(std::vector<VPL> &_vpls);

	LightNodeType *recursiveBuild(uint32_t nodeNum, int start, int end, 
		const VPL **vpls);

	LightNodeType *root;

	uint32_t nMaxNodes, nextFreeNode;
	std::vector<LightNodeType> m_nodes;

	// control weights between positon and direction, useful in surfaceLight
	Float m_ratio;
};

template <typename LightNodeType> 
void LightTree<LightNodeType>::build(std::vector<VPL> &_vpls) {
	nMaxNodes = _vpls.size() * 4;
	nextFreeNode = 1;
	m_nodes.resize(nMaxNodes);
	std::vector<const VPL*> vpls(_vpls.size(), NULL);
	for (size_t i = 0; i < _vpls.size(); i++) {
		vpls[i] = &_vpls[i];
	}
	root = recursiveBuild(0, 0, _vpls.size(), &vpls[0]);
}

template <typename LightNodeType>
LightNodeType* LightTree<LightNodeType>::recursiveBuild(uint32_t nodeNum, int start, int end, 
		const VPL **vpls) {
	if (start + 1 == end) {
		m_nodes[nodeNum].initLeaf();
		m_nodes[nodeNum].light = *vpls[start];
		return &m_nodes[nodeNum];
	}

	AABB6 bound;
	for (int i = start; i < end; i++) {

		bound.expandBy()
	}

}

MTS_NAMESPACE_END


#endif