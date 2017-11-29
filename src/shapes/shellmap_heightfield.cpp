#include <mitsuba/render/shape.h>
#include <mitsuba/render/skdtree.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sensor.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/timer.h>
#include "../volume/tetra.h"

MTS_NAMESPACE_BEGIN

namespace {
	/// Temporary storage for patch-ray intersections
	struct PatchIntersectionRecord {
		Point p;
		int x, y;
		int blockX, blockY;
    };

    struct ShellmapIntersectionRecord {
        Frame geoFrameBlock;
        Vector dpduBlock;
        Vector dpdvBlock;
        /*
        Point tex;
        Vector norm;
        TangentSpace tang;
        */
    };
};

class ShellmapHeightfield : public Shape {
public:
    ShellmapHeightfield(const Properties &props) : Shape(props) {
        m_objectToWorld = props.getTransform("toWorld", Transform());
        m_shellFilename = props.getString("shellFilename");

        // correspond to z=1 in texture space
        m_maxHeight = props.getFloat("maxHeight", -1.0);

        m_debug = false;
    }

    ShellmapHeightfield(Stream *stream, InstanceManager *manager) 
        : Shape(stream, manager) {
        m_objectToWorld = Transform(stream);
        m_shellFilename = stream->readString();
        m_maxHeight = stream->readFloat();
        m_block = static_cast<Shape *>(manager->getInstance(stream));
        m_meshBound = static_cast<Shape *>(manager->getInstance(stream));
        configure();
    }

    virtual ~ShellmapHeightfield() {
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Shape::serialize(stream, manager);
        m_objectToWorld.serialize(stream);
        stream->writeString(m_shellFilename);
        stream->writeFloat(m_maxHeight);
        manager->serialize(stream, m_block.get());
        manager->serialize(stream, m_meshBound.get());
    }

    AABB getAABB() const {
        return m_aabb;
    }

    Float getSurfaceArea() const {
        return 0;
    }

    size_t getPrimitiveCount() const {
        return 1;
    }

    size_t getEffectivePrimitiveCount() const {
        return m_block->getEffectivePrimitiveCount();
    }

    bool rayIntersect(const Ray &_ray, Float mint, Float maxt, Float &t, void *tmp) const {
        // world space
        Point enterPt, exitPt;
        Float nearT = mint, farT = maxt;
        if (!m_aabb.rayIntersect(_ray, nearT, farT, enterPt, exitPt))
            return false;

        // shellmap (object) space
        Ray ray;
        m_worldToObject(_ray, ray);
        
        if (m_debug) {
            Log(EInfo, "------------------------------------");
 		    Log(EInfo, "ray world: o = (%.3f, %.3f, %.3f), d = (%.3f, %.3f, %.3f)",
 		    	_ray.o.x, _ray.o.y, _ray.o.z, _ray.d.x, _ray.d.y, _ray.d.z);
 	    	Log(EInfo, "ray local: o = (%.3f, %.3f, %.3f), d = (%.3f, %.3f, %.3f)",
     			ray.o.x, ray.o.y, ray.o.z, ray.d.x, ray.d.y, ray.d.z);
            Log(EInfo, "mint = %.6f, maxt = %.6f", mint, maxt);
        }

        Intersection its;
        if (!m_kdtree->rayIntersect(ray, its))
            return false;

        t = 0;
        
        Point posShell(ray.o);
        posShell += ray.d * its.t;
        mint -= its.t;
        maxt -= its.t;
        t += its.t;

        if (m_debug) {
            Log(EInfo, "first hit boundary: %.6f", its.t);
        }

        // ray marching in shellmap space
        while (true) {
            ray.o = posShell;

            Float tFar;
            Point texNear, texFar;
            if (!m_shell.rayIntersectTetrahedron(ray, tFar, texNear, texFar))
                break;

            if (m_debug) {
 			    Log(EInfo, "===== tetra intersect =====");
 			    Log(EInfo, "texNear: (%.6f, %.6f, %.6f)", texNear.x, texNear.y, texNear.z);
 			    Log(EInfo, "texFar: (%.6f, %.6f, %.6f)", texFar.x, texFar.y, texFar.z);
                Log(EInfo, "tFar = %.6f", tFar);
            }

            clampTexPoint(texNear);
            clampTexPoint(texFar);

            Point blockNear, blockFar;
            blockNear = m_textureToData.transformAffine(texNear);
            blockFar = m_textureToData.transformAffine(texFar);

            if (m_debug) {
                Log(EInfo, "blockNear: (%.6f, %.6f, %.6f)", blockNear.x, blockNear.y, blockNear.z);
                Log(EInfo, "blockFar: (%.6f, %.6f, %.6f)", blockFar.x, blockFar.y, blockFar.z);
            }

            Ray rayBlock(blockNear, normalize(blockFar - blockNear), 0);
            Float tBlock;
            PatchIntersectionRecord tempBlock;
            Float mintBlock = std::max(Epsilon, mint);
            Float maxtBlock = std::min((blockFar - blockNear).length(), maxt);
            
            if (m_block->rayIntersect(rayBlock, mintBlock, maxtBlock, 
                tBlock, (void*)(&tempBlock))) {
                // intersection in block space
                Intersection itsBlock;
                m_block->fillIntersectionRecord(rayBlock, (void*)(&tempBlock), itsBlock);

                Float ratio = (itsBlock.p - blockNear).length() / (blockFar - blockNear).length();
                t += tFar * ratio;

                if (m_debug) {
                    Log(EInfo, "%.6f, %.6f", tFar, ratio);
                }

                Point itsPShell = m_worldToObject(_ray(t));
                Point tex;
                /*
                Vector norm;
                TangentSpace tang;
                */
                if (!m_shell.lookupPoint(itsPShell, tex)) {
                    break;
                }
                    
                if (tmp) {
                    ShellmapIntersectionRecord &temp = *((ShellmapIntersectionRecord*)tmp);             
                    temp.geoFrameBlock = itsBlock.geoFrame;
                    temp.dpduBlock = itsBlock.dpdu;
                    temp.dpdvBlock = itsBlock.dpdv;
                    /*
                    temp.tex = tex;
                    temp.norm = norm;
                    temp.tang = tang;
                    */
                }
                
                if (m_debug) {
                    Log(EInfo, "*** found intersection ***");
                    Log(EInfo, "ratio = %.6f, final dist = %.6f", ratio, t);
                }

                return true;
            }
            else {
                // move to the next tetrahedron
                posShell += ray.d * tFar;
                mint -= tFar;
                maxt -= tFar;
                t += tFar;
            }

            if (maxt < -Epsilon)
                break;
        }

        t = std::numeric_limits<Float>::infinity();
        return false;
    }

    bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
        Float t;
        return rayIntersect(ray, mint, maxt, t, NULL);
    }

    void fillIntersectionRecord(const Ray &ray, const void *tmp, Intersection &its) const {
        ShellmapIntersectionRecord &temp = *((ShellmapIntersectionRecord*)tmp);
        
        its.p = ray(its.t);

        Point posShell = m_worldToObject(its.p);
        Point tex;
        Vector norm;
        TangentSpace tang;
        bool flag = m_shell.lookupPoint(posShell, tex, norm, tang);

        if (!flag) {
            Log(EInfo, "p = (%.6f, %.6f, %.6f), dist = %.6f", its.p.x, its.p.y, its.p.z, its.t);
            Log(EError, "Not inside tetrahedra!");
            return;
        }

        /*
        Point &tex = temp.tex;
        Vector &norm = temp.norm;
        TangentSpace &tang = temp.tang;
        */

        its.uv = Point2(tex.x, tex.y);

        Vector dpduTex = m_dataToTexture(temp.dpduBlock);
        Vector dpdvTex = m_dataToTexture(temp.dpdvBlock);
        Normal normTex = normalize(m_dataToTexture(temp.geoFrameBlock.n));
       
        Vector dpduShell = dpduTex.x * tang.dpdu + dpduTex.y * tang.dpdv + dpduTex.z * norm;
        Vector dpdvShell = dpdvTex.x * tang.dpdu + dpdvTex.y * tang.dpdv + dpdvTex.z * norm;
        Normal normShell = normTex.x * tang.dpdu + normTex.y * tang.dpdv + normTex.z * norm;

        Vector normWorld = normalize(m_objectToWorld(normShell));

        its.dpdu = m_objectToWorld(dpduShell);
        its.dpdv = m_objectToWorld(dpdvShell);

        /*
        its.geoFrame.n = normWorld;
        its.geoFrame.s = normalize(its.dpdu);
        its.geoFrame.t = cross(its.geoFrame.n, its.geoFrame.s);
        */

        its.geoFrame.s = normalize(its.dpdu);
        its.geoFrame.t = normalize(its.dpdv - dot(its.dpdv, its.geoFrame.s) * its.geoFrame.s);
        its.geoFrame.n = cross(its.geoFrame.s, its.geoFrame.t);

        its.shFrame.n = its.geoFrame.n;

        its.shape = this;
        its.hasUVPartials = false;
        its.instance = NULL;
        its.time = ray.time;
    }

    //void getNormalDerivative(const Intersection &its, Vector &dndu, Vector &dndv, 
    //    bool shadingFrame) const {
    //}

    void configure() {
        if (m_block.get() == NULL) 
            Log(EError, "No embedded heightfield specified!");
        if (m_meshBound.get() == NULL)
            Log(EError, "No mesh boundary specified!");

        m_worldToObject = m_objectToWorld.inverse();
        
        AABB blockAABB = m_block->getAABB();
        Float maxz = m_maxHeight;
        /*
        if (blockAABB.max.z - blockAABB.min.z > Epsilon) {
            maxz = blockAABB.max.z + blockAABB.getExtents().z * 0.01;
        }
        else {
            if (m_maxHeight < 0)
                Log(EError, "Singular texture mapping in z!");
            maxz = m_maxHeight;
        }
        */
        blockAABB.min.z = 0;
        blockAABB.max.z = maxz;
        
        m_textureToData = Transform::translate(Vector(blockAABB.min)) * 
            Transform::scale(blockAABB.getExtents());

        m_dataToTexture = m_textureToData.inverse();

        fs::path resolved = Thread::getThread()->getFileResolver()->resolve(m_shellFilename);
        if (!m_shell.load(resolved.string().c_str()))
            Log(EError, "Failed to load the shell file!");
        else
            Log(EInfo, "Shell mesh loaded: %u tetrahedra, tree depth: %u",
                m_shell.getTetrahedronCount(), m_shell.getTreeDepth());

        m_aabb.reset();
        for (int i = 0; i < 8; ++i)
            m_aabb.expandBy(m_objectToWorld(m_meshBound->getAABB().getCorner(i)));

        // build kd-tree of the mesh boundary
        m_kdtree = new ShapeKDTree();
        if (!m_kdtree->isBuilt()) {
            addShape(m_meshBound);
            m_kdtree->build();
        }
        Log(EInfo, "Shellmap heightfield configuration finished");
        Log(EInfo, "%s", toString().c_str());
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Shape)) && name == "baseHeightfield") {
            Assert(m_block == NULL);
            m_block = static_cast<Shape*>(child);
        }
        else if (child->getClass()->derivesFrom(MTS_CLASS(Shape)) && name == "meshBound") {
            Assert(m_meshBound == NULL);
            m_meshBound = static_cast<Shape*>(child);
        }
        else {
            Shape::addChild(name, child);
        }
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "ShellmapHeightfield[" << endl
            << "  objectToWorld = " << indent(m_objectToWorld.toString()) << "," << endl
            << "  textureToData = " << indent(m_textureToData.toString()) << "," << endl
            << "  aabb = " << indent(getAABB().toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()

private:
    void clampTexPoint(Point &p) const {
        if (p.z < -Epsilon || p.z > 1.0f + Epsilon) {
            //Log(EError, "bad z = %.6f", p.z);
            p.z = math::clamp(p.z, 0.0, 1.0);
        }
        p.x -= math::floorToInt(p.x);
        p.y -= math::floorToInt(p.y);
    }

    void addShape(Shape *shape) {
        if (shape->isCompound()) {
            int index = 0;
            do {
                ref<Shape> element = shape->getElement(index++);
                if (element == NULL)
                    break;
                addShape(element);
            } while (true);
        }
        else {
            m_kdtree->addShape(shape);
        }
    }

protected:
    std::string m_shellFilename;
    ref<Shape> m_block;
    ref<Shape> m_meshBound;
    ref<ShapeKDTree> m_kdtree;
    TetrahedronMesh m_shell;
    Float m_maxHeight;
    Transform m_worldToObject, m_objectToWorld;
    Transform m_textureToData, m_dataToTexture;
    AABB m_aabb;

    bool m_debug;
};

MTS_IMPLEMENT_CLASS_S(ShellmapHeightfield, false, Shape)
MTS_EXPORT_PLUGIN(ShellmapHeightfield, "Shellmapped height field intersection shape");
MTS_NAMESPACE_END
