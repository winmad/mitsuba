#ifdef __TETRAHEDRON_H
#error "badness"
#endif

#ifndef __TETRAHEDRON2_H
#define __TETRAHEDRON2_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/render/skdtree.h>
#include <mitsuba/core/tls.h>
#include <vector>
#include <map>


MTS_NAMESPACE_BEGIN


struct Tetrahedron
{
    uint32_t idx[4];

    inline AABB getAABB(const Point *positions) const
    {
        AABB result(positions[idx[0]]);
        result.expandBy(positions[idx[1]]);
        result.expandBy(positions[idx[2]]);
        result.expandBy(positions[idx[3]]);
        return result;
    }

    inline bool inside(const Point *positions, const Point &p, Point4 &bb) const
    {
        Vector3 a(positions[idx[0]] - p);
        Vector3 b(positions[idx[1]] - p);
        Vector3 c(positions[idx[2]] - p);
        Vector3 d(positions[idx[3]] - p);

        Float iv = 1.0f/std::abs(dot(a - d, cross(b - d, c - d)));
        Float v0 = std::abs(dot(b, cross(c, d)));
        Float v1 = std::abs(dot(a, cross(c, d)));
        Float v2 = std::abs(dot(a, cross(b, d)));
        Float v3 = std::abs(dot(a, cross(b, c)));

        bb.x = v0*iv;
        bb.y = v1*iv;
        bb.z = v2*iv;
        bb.w = v3*iv;

        return bb.x + bb.y + bb.z + bb.w < 1.0f + Epsilon;
    }

    bool rayIntersect(const Point *positions, const Ray &r, Float &minT, uint32_t &minFace, Float &maxT, uint32_t &maxFace) const
    {
        //Ray r(_r);
        //r.mint = Epsilon;
        //r.maxt = std::numeric_limits<Float>::infinity();

        const int loc[] = {0, 1, 2, 0, 2, 3, 0, 3, 1, 1, 3, 2};
        Triangle tri[4];
        for ( uint32_t i = 0; i < 4; ++i )
        {
            tri[i].idx[0] = idx[loc[3*i]];
            tri[i].idx[1] = idx[loc[3*i + 1]];
            tri[i].idx[2] = idx[loc[3*i + 2]];
        }

        //tri[0].idx[0] = idx[0]; tri[0].idx[1] = idx[1]; tri[0].idx[2] = idx[2];
        //tri[1].idx[0] = idx[0]; tri[1].idx[1] = idx[2]; tri[1].idx[2] = idx[3];
        //tri[2].idx[0] = idx[0]; tri[2].idx[1] = idx[3]; tri[2].idx[2] = idx[1];
        //tri[3].idx[0] = idx[1]; tri[3].idx[1] = idx[3]; tri[3].idx[2] = idx[2];

        minT = std::numeric_limits<Float>::infinity();
        maxT = -std::numeric_limits<Float>::infinity();

        for ( uint32_t i = 0; i < 4; ++i )
        {
            Float u, v, t;
            if ( tri[i].rayIntersect(positions, r, u, v, t) )
            {
                Vector e1 = positions[tri[i].idx[1]] - positions[tri[i].idx[0]];
                Vector e2 = positions[tri[i].idx[2]] - positions[tri[i].idx[0]];
                Vector n = normalize(cross(e1, e2));

                Float v = dot(n, r.d);
                if ( v < 0.0f )
                {
                    minT = t; minFace = i;
                }
                else //if ( v > Epsilon )
                {
                    maxT = t; maxFace = i;
                }
            }
        }

        if ( maxT - minT > -Epsilon )
            return true;
        else
        {
            //std::cout << minT << ' ' << maxT << std::endl;
            return false;
        }
    }
};


class TetrahedronMesh
{
public:
    TetrahedronMesh()
    {
        m_vertexCount = m_tetrahedronCount = 0;
        m_vtxPosition = m_vtxTexcoord = NULL;
        m_vtxNormal = NULL; m_vtxTangent = NULL;
        m_tetra = NULL;

        m_list = NULL; m_tree = NULL;
        m_treeSize = 0;

        m_list2 = NULL; m_tree2 = NULL;
        m_treeSize2 = 0;
    }

    ~TetrahedronMesh()
    {
        if ( m_vtxPosition ) delete[] m_vtxPosition;
        if ( m_vtxTexcoord ) delete[] m_vtxTexcoord;
        if ( m_vtxNormal ) delete[] m_vtxNormal;
        if ( m_vtxTangent ) delete[] m_vtxTangent;
        if ( m_tetra ) delete[] m_tetra;

        if ( m_list ) delete[] m_list;
        if ( m_tree ) delete[] m_tree;

        if ( m_list2 ) delete[] m_list2;
        if ( m_tree2 ) delete[] m_tree2;
    }

    bool load(const char *fname)
    {
        #define VERIFY_VALUE(x, y)	{ if ( (x) != (y) ) return false; }

        FILE *fin = fopen(fname, "rt");
        if ( fin )
        {
            VERIFY_VALUE(fscanf(fin, "%u %u", &m_vertexCount, &m_tetrahedronCount), 2);

            m_vtxPosition = new Point[m_vertexCount];
            m_vtxTexcoord = new Point[m_vertexCount];
            m_vtxNormal = new Vector[m_vertexCount];
            m_vtxTangent = new TangentSpace[m_vertexCount];
            m_tetra = new Tetrahedron[m_tetrahedronCount];

            for ( uint32_t i = 0; i < m_vertexCount; ++i )
            {
#ifdef DOUBLE_PRECISION
				VERIFY_VALUE(fscanf(fin, "%lf %lf %lf", &m_vtxPosition[i].x, &m_vtxPosition[i].y, &m_vtxPosition[i].z), 3);
				VERIFY_VALUE(fscanf(fin, "%lf %lf %lf", &m_vtxTexcoord[i].x, &m_vtxTexcoord[i].y, &m_vtxTexcoord[i].z), 3);
				VERIFY_VALUE(fscanf(fin, "%lf %lf %lf", &m_vtxNormal[i].x, &m_vtxNormal[i].y, &m_vtxNormal[i].z), 3);
				VERIFY_VALUE(fscanf(fin, "%lf %lf %lf", &m_vtxTangent[i].dpdu.x, &m_vtxTangent[i].dpdu.y, &m_vtxTangent[i].dpdu.z), 3);
				VERIFY_VALUE(fscanf(fin, "%lf %lf %lf", &m_vtxTangent[i].dpdv.x, &m_vtxTangent[i].dpdv.y, &m_vtxTangent[i].dpdv.z), 3);
#else
				VERIFY_VALUE(fscanf(fin, "%f %f %f", &m_vtxPosition[i].x, &m_vtxPosition[i].y, &m_vtxPosition[i].z), 3);
				VERIFY_VALUE(fscanf(fin, "%f %f %f", &m_vtxTexcoord[i].x, &m_vtxTexcoord[i].y, &m_vtxTexcoord[i].z), 3);
				VERIFY_VALUE(fscanf(fin, "%f %f %f", &m_vtxNormal[i].x, &m_vtxNormal[i].y, &m_vtxNormal[i].z), 3);
				VERIFY_VALUE(fscanf(fin, "%f %f %f", &m_vtxTangent[i].dpdu.x, &m_vtxTangent[i].dpdu.y, &m_vtxTangent[i].dpdu.z), 3);
				VERIFY_VALUE(fscanf(fin, "%f %f %f", &m_vtxTangent[i].dpdv.x, &m_vtxTangent[i].dpdv.y, &m_vtxTangent[i].dpdv.z), 3);
#endif

// 				Float cosineTerm = dot(m_vtxTangent[i].dpdu, m_vtxTangent[i].dpdv);
// 				if (std::abs(cosineTerm) > 1e-4) {
// 					printf("%.6f, Not an orthogonal frame (input)\n", cosineTerm);
// 				}
            }
            for ( uint32_t i = 0; i < m_tetrahedronCount; ++i )
                VERIFY_VALUE(fscanf(fin, "%u %u %u %u", &m_tetra[i].idx[0], &m_tetra[i].idx[1], &m_tetra[i].idx[2], &m_tetra[i].idx[3]), 4);

            fclose(fin);

            return true;
        }
        else
            return false;

        #undef VERIFY_VALUE
    }

    void configure(Transform toWorld = Transform(), bool buildInverseTree = false)
    {
        // transform the mesh into the world space
        for ( uint32_t i = 0; i < m_vertexCount; ++i )
        {
            m_vtxPosition[i] = toWorld.transformAffine(m_vtxPosition[i]);
            m_vtxNormal[i] = toWorld(m_vtxNormal[i]);
            m_vtxTangent[i].dpdu = toWorld(m_vtxTangent[i].dpdu);
            m_vtxTangent[i].dpdv = toWorld(m_vtxTangent[i].dpdv);
        }

        // compute AABB
        m_aabb.reset();
        for ( uint32_t i = 0; i < m_tetrahedronCount; ++i )
            m_aabb.expandBy(m_tetra[i].getAABB(m_vtxPosition));

        // fix vertex ordering for each face
        for ( uint32_t i = 0; i < m_tetrahedronCount; ++i )
        {
            Tetrahedron &t = m_tetra[i];
            Vector e1 = normalize(m_vtxPosition[t.idx[2]] - m_vtxPosition[t.idx[1]]);
            Vector e2 = normalize(m_vtxPosition[t.idx[3]] - m_vtxPosition[t.idx[1]]);
            Vector e = cross(e1, e2);
            if ( dot(e, Vector(m_vtxPosition[t.idx[0]])) < dot(e, Vector(m_vtxPosition[t.idx[1]])) )
                std::swap(t.idx[2], t.idx[3]);
        }

        // build BVH
        {
            m_list = new uint32_t[m_tetrahedronCount];
            m_tree = new _node_type[2*m_tetrahedronCount];

            AABB *abl = new AABB[m_tetrahedronCount], aabb;
            for ( uint32_t i = 0; i < m_tetrahedronCount; ++i )
            {
                m_list[i] = i;
                abl[i] = m_tetra[i].getAABB(m_vtxPosition);
                aabb.expandBy(abl[i]);
            }
            std::random_shuffle(m_list, m_list + m_tetrahedronCount);

            _node_type &root = m_tree[1];
            root.L = 0; root.R = m_tetrahedronCount;
            root.Lch = root.Rch = -1;
            root.aabb = aabb;
            m_treeSize = 1;
            m_treeDepth = build(0, 1, abl, m_treeSize, m_tree, m_list);
            SAssert( m_treeDepth <= MAX_TREE_DEPTH );

            delete[] abl;
        }

        if ( buildInverseTree )
        {
            m_list2 = new uint32_t[m_tetrahedronCount];
            m_tree2 = new _node_type[2*m_tetrahedronCount];

            AABB *abl = new AABB[m_tetrahedronCount], aabb;
            for ( uint32_t i = 0; i < m_tetrahedronCount; ++i )
            {
                m_list2[i] = i;
                abl[i] = m_tetra[i].getAABB(m_vtxTexcoord);
                aabb.expandBy(abl[i]);
            }
            std::random_shuffle(m_list2, m_list2 + m_tetrahedronCount);

            _node_type &root2 = m_tree2[1];
            root2.L = 0; root2.R = m_tetrahedronCount;
            root2.Lch = root2.Rch = -1;
            root2.aabb = aabb;
            m_treeSize2 = 1;
            m_treeDepth2 = build(0, 1, abl, m_treeSize2, m_tree2, m_list2);
            SAssert( m_treeDepth2 <= MAX_TREE_DEPTH );

            delete[] abl;
        }
        else
        {
            m_list2 = NULL;
            m_tree2 = NULL;
            m_treeDepth2 = m_treeSize2 = 0;
        }

        // build indices
        _map_type map;
        for ( uint32_t i = 0; i < m_tetrahedronCount; ++i )
        {
            const Tetrahedron &t = m_tetra[i];

            _face_type f[4];
            f[0] = _face_type(t.idx[0], t.idx[1], t.idx[2]);
            f[1] = _face_type(t.idx[0], t.idx[2], t.idx[3]);
            f[2] = _face_type(t.idx[0], t.idx[3], t.idx[1]);
            f[3] = _face_type(t.idx[1], t.idx[3], t.idx[2]);

            for ( uint32_t j = 0; j < 4; ++j )
                map.insert(_map_type::value_type(f[j], i));
        }

        m_link.resize(4*m_tetrahedronCount);

        for ( uint32_t i = 0; i < m_tetrahedronCount; ++i )
        {
            Tetrahedron &t = m_tetra[i];

            _face_type f[4];
            f[0] = _face_type(t.idx[0], t.idx[1], t.idx[2]);
            f[1] = _face_type(t.idx[0], t.idx[2], t.idx[3]);
            f[2] = _face_type(t.idx[0], t.idx[3], t.idx[1]);
            f[3] = _face_type(t.idx[1], t.idx[3], t.idx[2]);

            _map_type::const_iterator it;
            for ( uint32_t j = 0; j < 4; ++j )
            {
                _face_type bf(f[j]); bf.reverse();
                if ( (it = map.find(bf)) != map.end() )
                    m_link[4*i + j] = static_cast<int>(it->second);
                else
                {
                    m_link[4*i + j] = -1;

                    m_boundaryFace.push_back(f[j]);
                    m_boundaryID.push_back(4*i + j);
                }
            }
        }

        ref<TriMesh> m_bmesh = new TriMesh("Test", m_boundaryFace.size(), 3*m_boundaryFace.size());
        Point *vtx = m_bmesh->getVertexPositions();
        Triangle *tri = m_bmesh->getTriangles();

        uint32_t tot = 0;
        for ( std::vector<_face_type>::const_iterator it = m_boundaryFace.begin(); it != m_boundaryFace.end(); ++it )
        {
            vtx[3*tot    ] = m_vtxPosition[it->idx[0]];
            vtx[3*tot + 1] = m_vtxPosition[it->idx[1]];
            vtx[3*tot + 2] = m_vtxPosition[it->idx[2]];

            tri[tot].idx[0] = 3*tot;
            tri[tot].idx[1] = 3*tot + 1;
            tri[tot].idx[2] = 3*tot + 2;

            ++tot;
        }
        m_bmesh->configure();
        m_surfaceArea = m_bmesh->getSurfaceArea();

        m_btree = new ShapeKDTree();
        m_btree->addShape(m_bmesh);
        m_btree->build();
    }

    inline AABB getAABB() const
    {
        return m_aabb;
    }

    inline Float getSurfaceArea() const
    {
        return m_surfaceArea;
    }

    bool lookupPoint_BruteForce(const Point &p, Point &tex) const
    {
        Point4 bb;
        for ( uint32_t i = 0; i < m_tetrahedronCount; ++i )
            if ( m_tetra[i].inside(m_vtxPosition, p, bb) )
            {
                tex = m_vtxTexcoord[m_tetra[i].idx[0]]*bb.x
                    + m_vtxTexcoord[m_tetra[i].idx[1]]*bb.y
                    + m_vtxTexcoord[m_tetra[i].idx[2]]*bb.z
                    + m_vtxTexcoord[m_tetra[i].idx[3]]*bb.w;
                return true;
            }
        return false;
    }

    bool lookupPoint(const Point &p, Point &tex) const
    {
        uint32_t id;
        Point4 bb;

        if ( !lookup(p, id, bb) ) return false;
        const uint32_t *idx = m_tetra[id].idx;
        tex = m_vtxTexcoord[idx[0]]*bb.x + m_vtxTexcoord[idx[1]]*bb.y
            + m_vtxTexcoord[idx[2]]*bb.z + m_vtxTexcoord[idx[3]]*bb.w;

        return true;
    }

    bool lookupPoint2(const Point &tex, Point &p) const
    {
        uint32_t id;
        Point4 bb;

        if ( !lookup2(tex, id, bb) ) return false;
        const uint32_t *idx = m_tetra[id].idx;
        p = m_vtxPosition[idx[0]]*bb.x + m_vtxPosition[idx[1]]*bb.y
          + m_vtxPosition[idx[2]]*bb.z + m_vtxPosition[idx[3]]*bb.w;
        return true;
    }

    bool lookupPoint(const Point &p, uint32_t &id) const
    {
        return lookup(p, id);
    }

    bool lookupPoint(const Point &p, Point &tex, Vector &normal, TangentSpace &tangent) const
    {
        uint32_t id;
        Point4 bb;
        
        //if ( !lookup(p, 1, id, bb) ) return false;
        if ( !lookup(p, id, bb) ) return false;
        const uint32_t *idx = m_tetra[id].idx;
        tex = m_vtxTexcoord[idx[0]]*bb.x + m_vtxTexcoord[idx[1]]*bb.y
            + m_vtxTexcoord[idx[2]]*bb.z + m_vtxTexcoord[idx[3]]*bb.w;
        normal = m_vtxNormal[idx[0]]*bb.x + m_vtxNormal[idx[1]]*bb.y
               + m_vtxNormal[idx[2]]*bb.z + m_vtxNormal[idx[3]]*bb.w;
        tangent.dpdu = m_vtxTangent[idx[0]].dpdu*bb.x + m_vtxTangent[idx[1]].dpdu*bb.y
                     + m_vtxTangent[idx[2]].dpdu*bb.z + m_vtxTangent[idx[3]].dpdu*bb.w;
        tangent.dpdv = m_vtxTangent[idx[0]].dpdv*bb.x + m_vtxTangent[idx[1]].dpdv*bb.y
                     + m_vtxTangent[idx[2]].dpdv*bb.z + m_vtxTangent[idx[3]].dpdv*bb.w;

        return true;
    }

    bool lookupPointGivenId(const Point &p, uint32_t id, Point *tex = NULL,
        Vector *normal = NULL, TangentSpace *tangent = NULL) const {
        Point4 bb;

        if (!m_tetra[id].inside(m_vtxPosition, p, bb)) 
            return false;
        const uint32_t *idx = m_tetra[id].idx;

        if (tex != NULL) {
            *tex = m_vtxTexcoord[idx[0]]*bb.x + m_vtxTexcoord[idx[1]]*bb.y
                + m_vtxTexcoord[idx[2]]*bb.z + m_vtxTexcoord[idx[3]]*bb.w;
        }

        if (normal != NULL) {
            *normal = m_vtxNormal[idx[0]]*bb.x + m_vtxNormal[idx[1]]*bb.y
                + m_vtxNormal[idx[2]]*bb.z + m_vtxNormal[idx[3]]*bb.w;
        }
            
        if (tangent != NULL) {
            tangent->dpdu = m_vtxTangent[idx[0]].dpdu*bb.x + m_vtxTangent[idx[1]].dpdu*bb.y
                + m_vtxTangent[idx[2]].dpdu*bb.z + m_vtxTangent[idx[3]].dpdu*bb.w;
            tangent->dpdv = m_vtxTangent[idx[0]].dpdv*bb.x + m_vtxTangent[idx[1]].dpdv*bb.y
                + m_vtxTangent[idx[2]].dpdv*bb.z + m_vtxTangent[idx[3]].dpdv*bb.w;
        }
        return true;
    }

    inline uint32_t getTetrahedronCount() const { return m_tetrahedronCount; }

    inline uint32_t getTreeSize() const         { return m_treeSize; }
    inline uint32_t getTreeDepth() const        { return m_treeDepth; }
    inline uint32_t getTreeSize2() const        { return m_treeSize2; }
    inline uint32_t getTreeDepth2() const       { return m_treeDepth2; }

protected:
    struct _node_type
    {
        _node_type()
        {
            L = R = 0;
            Lch = Rch = -1;
        }

        _node_type(const _node_type& in)
        {
            L = in.L; R = in.R;
            Lch = in.Lch; Rch = in.Rch;
            aabb = in.aabb;
        }

        uint32_t L, R;
        int Lch, Rch;
        AABB aabb;
    };

    static const uint32_t INVALID_INDEX = 0xffffffff;
    static const uint32_t MAX_TREE_DEPTH = 100;

    static uint32_t build(uint32_t dep, uint32_t r, const AABB *abl,
        uint32_t &treeSize, _node_type *tree, uint32_t *list)
    {
        if ( dep == MAX_TREE_DEPTH /*|| tree[r].R - tree[r].L <= 5*/ )
            return 1;

        int axis = tree[r].aabb.getLargestAxis();
        Float pos = 0.5f*(tree[r].aabb.min[axis] + tree[r].aabb.max[axis]);

        int L = tree[r].L, R = tree[r].R;
        while ( L < R )
        {
            while ( L < R && pos - abl[list[L]].getCenter()[axis] > -Epsilon ) ++L;
            while ( L < R && abl[list[R - 1]].getCenter()[axis] - pos > -Epsilon ) --R;
            if ( L < R ) std::swap(list[L], list[R - 1]);
        }

        if ( L == static_cast<int>(tree[r].L) || L == static_cast<int>(tree[r].R) )
            return 1;
        
        _node_type nodeL, nodeR;
        int i;

        nodeL.L = tree[r].L; nodeL.R = L;
        nodeL.Lch = nodeL.Rch = -1;
        for ( i = tree[r].L; i < L; ++i )
            nodeL.aabb.expandBy(abl[list[i]]);
        nodeR.L = L; nodeR.R = tree[r].R;
        nodeR.Lch = nodeR.Rch = -1;
        for ( i = static_cast<int>(L); i < static_cast<int>(tree[r].R); ++i )
            nodeR.aabb.expandBy(abl[list[i]]);

        tree[++treeSize] = nodeL;
        tree[++treeSize] = nodeR;
        tree[r].Lch = static_cast<int>(treeSize) - 1;
        tree[r].Rch = static_cast<int>(treeSize);
        return std::max(build(dep + 1, tree[r].Lch, abl, treeSize, tree, list),
            build(dep + 1, tree[r].Rch, abl, treeSize, tree, list)) + 1;
    }

    static uint32_t lookup(const Point &p, uint32_t r, uint32_t &id, Point4 &bb,
        const _node_type *tree, const uint32_t *list, const Tetrahedron *tetra, const Point *vtxPosition,
        PrimitiveThreadLocal<uint32_t> &cache_tetra)
    {
        if ( !tree[r].aabb.contains(p) ) return 0;

        uint32_t j;
        if ( tree[r].Lch < 0 )
        {
            for ( uint32_t i = tree[r].L; i < tree[r].R; ++i )
            {
                j = list[i];
                if ( tetra[j].inside(vtxPosition, p, bb) )
                {
                    id = j;
#ifndef TETRAHEDRON_MESH_NO_CACHE
                    cache_tetra.get() = j;
#endif
                    return 1;
                }
            }
        }
        else
        {
            uint32_t ret;

            j = tree[r].Lch;
            ret = lookup(p, j, id, bb, tree, list, tetra, vtxPosition, cache_tetra);
            if ( ret ) return ret + 1;

            j = tree[r].Rch;
            ret = lookup(p, j, id, bb, tree, list, tetra, vtxPosition, cache_tetra);
            if ( ret ) return ret + 1;
        }

        return 0;
    }

    bool lookup(const Point &p, uint32_t &id, Point4 &bb) const
    {
#ifndef TETRAHEDRON_MESH_NO_CACHE
        uint32_t val = m_cache_tetra.get();
        if ( val <= m_treeSize && m_tetra[val].inside(m_vtxPosition, p, bb) )
        {
            id = val; return true;
        }
#endif
        return lookup(p, 1, id, bb, m_tree, m_list, m_tetra, m_vtxPosition, m_cache_tetra);
    }

    bool lookup(const Point &p, uint32_t &id) const
    {
        Point4 bb;
        return lookup(p, id, bb);
    }

    bool lookup2(const Point &tex, uint32_t &id, Point4 &bb) const
    {
        if ( m_tree2 == NULL )
        {
            std::cout << "Tetrahedron Mesh: backward mapping not available." << std::endl;
            return false;
        }
#ifndef TETRAHEDRON_MESH_NO_CACHE
        uint32_t val = m_cache_tetra2.get();
        if ( val <= m_treeSize2 && m_tetra[val].inside(m_vtxTexcoord, tex, bb) )
        {
            id = val; return true;
        }
#endif
        return lookup(tex, 1, id, bb, m_tree2, m_list2, m_tetra, m_vtxTexcoord, m_cache_tetra2);
    }

    bool lookup2(const Point &tex, uint32_t &id) const
    {
        Point4 bb;
        return lookup(tex, id, bb);
    }

public:
    AABB m_aabb;
    Float m_surfaceArea;

    uint32_t m_vertexCount, m_tetrahedronCount;
    Point *m_vtxPosition, *m_vtxTexcoord;
    Vector *m_vtxNormal;
    TangentSpace *m_vtxTangent;
    Tetrahedron *m_tetra;

protected:
    // BVH for point -> tex coord lookups
    uint32_t *m_list;
    uint32_t m_treeSize, m_treeDepth;
    _node_type *m_tree;

    // BVH for tex coord -> point lookups
    uint32_t *m_list2;
    uint32_t m_treeSize2, m_treeDepth2;
    _node_type *m_tree2;

    // per-thread cache
#ifdef TETRAHEDRON_MESH_NO_CACHE
    #pragma message("Tetrahedron mesh compiled without per-thread cache")
#endif
    mutable PrimitiveThreadLocal<uint32_t> m_cache_tetra;
    mutable PrimitiveThreadLocal<uint32_t> m_cache_tetra2;

    struct _face_type
    {
        _face_type()
        {
            memset(idx, 0, sizeof(idx));
        }

        _face_type(uint32_t v0, uint32_t v1, uint32_t v2)
        {
            uint32_t v[3];
            v[0] = v0; v[1] = v1; v[2]= v2;

            int i = 0;
            for ( int j = 1; j < 3; ++j )
                if ( v[j] < v[i] ) i = j;

            for ( int j = 0; j < 3; ++j )
                idx[j] = v[(i + j) % 3];
        }

        _face_type(const _face_type &q)
        {
            memcpy(idx, q.idx, 3*sizeof(uint32_t));
        }

        bool operator< (const _face_type &q) const
        {
            for ( int i = 0; i < 3; ++i )
                if ( idx[i] != q.idx[i] ) return idx[i] < q.idx[i];
            return false;
        }

        void reverse()
        {
            std::swap(idx[1], idx[2]);
        }

        uint32_t idx[3];
    };

    typedef std::map<_face_type, uint32_t> _map_type;

public:
    std::vector<int> m_link;
    std::vector<_face_type> m_boundaryFace;
    std::vector<uint32_t> m_boundaryID;

    ref<ShapeKDTree> m_btree;
};


MTS_NAMESPACE_END


#endif
