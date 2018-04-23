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
#include "../volume/tetra2.h"

//#define SHELLMAP_HEIGHTFIELD_DEBUG

MTS_NAMESPACE_BEGIN

namespace {
	/// Temporary storage for patch-ray intersections
	struct PatchIntersectionRecord {
		Point p;
		int x, y;
		int blockX, blockY;
    };

    struct ShellmapIntersectionRecord {
        Normal normBlock;
        Vector dpduBlock;
        Vector dpdvBlock;
        Point4 bb;
        int tetraId;
    };
};

class ShellmapHeightfield : public Shape {
public:
    ShellmapHeightfield(const Properties &props) : Shape(props), m_ready(false) {
        m_objectToWorld = props.getTransform("toWorld", Transform());
        m_shellFilename = props.getString("shellFilename");

        // correspond to z=1 in texture space
        m_maxHeight = props.getFloat("maxHeight", -1.0);

		m_useMacroDeform = props.getBoolean("useMacroDeform", false);
    }

    ShellmapHeightfield(Stream *stream, InstanceManager *manager) 
        : Shape(stream, manager), m_ready(false) {
        m_objectToWorld = Transform(stream);
        m_shellFilename = stream->readString();
        m_maxHeight = stream->readFloat();
		m_useMacroDeform = stream->readBool();
        m_block = static_cast<Shape *>(manager->getInstance(stream));
        configure();
    }

    virtual ~ShellmapHeightfield() {
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Shape::serialize(stream, manager);
        m_objectToWorld.serialize(stream);
        stream->writeString(m_shellFilename);
        stream->writeFloat(m_maxHeight);
		stream->writeBool(m_useMacroDeform);
        manager->serialize(stream, m_block.get());
    }

    AABB getAABB() const {
        return m_aabb;
    }

    Float getSurfaceArea() const {
        return m_shell.getSurfaceArea();
    }

    size_t getPrimitiveCount() const {
        return 1;
    }

    size_t getEffectivePrimitiveCount() const {
        return m_block->getEffectivePrimitiveCount();
    }

    bool rayIntersect(const Ray &ray, Float mint, Float maxt, Float &t, void *tmp) const {
        Float _mint, _maxt;
        Ray _ray;

        uint32_t id0, f0, f1;
        Float t0, t1;
        if (m_shell.lookupPoint(ray(mint), id0)) {
            _mint = mint; 
            _maxt = maxt;
        }
        else {
            Intersection its;
            _ray = Ray(ray, mint, std::numeric_limits<Float>::infinity());
            if (!m_shell.m_btree->rayIntersect(_ray, its) || its.t > maxt)
                return false;
            
            _mint = its.t;
            if (!m_shell.lookupPoint(ray(_mint + Epsilon), id0))
                return false;

            _maxt = maxt;
        }

        _ray = Ray(ray, -std::numeric_limits<Float>::infinity(), std::numeric_limits<Float>::infinity());
        if (!m_shell.m_tetra[id0].rayIntersect(m_shell.m_vtxPosition, _ray, t0, f0, t1, f1)) {
            return false;
        }

        Float tTemp = _mint;
        int id = static_cast<int>(id0);
        Float minT = t0, maxT = t1;
        uint32_t minF = f0, maxF = f1;

#ifdef SHELLMAP_HEIGHTFIELD_DEBUG
        Log(EInfo, "==========================");
        Log(EInfo, "ray world: o = (%.6f, %.6f, %.6f), d = (%.6f, %.6f, %.6f)",
 			ray.o.x, ray.o.y, ray.o.z, ray.d.x, ray.d.y, ray.d.z);
        Log(EInfo, "init mint = %.8f, maxt = %.8f", mint, maxt);
        Log(EInfo, "intersect tetra %d: minT = %.8f, maxT = %.8f", id, minT, maxT);
#endif

		int cnt = 0;

        while (true) {
            if (tTemp >= _maxt)
                break;

			if (++cnt > m_tetrahedronCount) {
				Log(EInfo, "Loop forever when intersecting...");
				break;
			}

#ifdef SHELLMAP_HEIGHTFIELD_DEBUG
            Log(EInfo, "-----------------------");
            Log(EInfo, "tetra id = %d, t = %.8f", id, tTemp);
            Log(EInfo, "intersect tetra %d: minT = %.8f, maxT = %.8f", id, minT, maxT);
#endif

            minT = std::max(minT, _mint);
            maxT = std::min(maxT, _maxt);
            
			if (minT < maxT) {
				Point texNear, texFar;
				m_shell.lookupPointGivenId(_ray(minT + Epsilon), id, &texNear);
				m_shell.lookupPointGivenId(_ray(maxT - Epsilon), id, &texFar);
				clampTexPoint(texNear);
				clampTexPoint(texFar);

				Point blockNear, blockFar;
				blockNear = m_textureToData.transformAffine(texNear);
				blockFar = m_textureToData.transformAffine(texFar);
				Float blockLength = (blockFar - blockNear).length();
				Ray rayBlock = Ray(blockNear, normalize(blockFar - blockNear),
					Epsilon, blockLength, 0);

#ifdef SHELLMAP_HEIGHTFIELD_DEBUG
				Log(EInfo, "World: near = (%.6f, %.6f, %.6f), far = (%.6f, %.6f, %.6f)",
					_ray(minT).x, _ray(minT).y, _ray(minT).z,
					_ray(maxT).x, _ray(maxT).y, _ray(maxT).z);
				Log(EInfo, "Tex: near = (%.6f, %.6f, %.6f), far = (%.6f, %.6f, %.6f)",
					texNear.x, texNear.y, texNear.z,
					texFar.x, texFar.y, texFar.z);
				Log(EInfo, "intersection length in block = %.8f", blockLength);
				Log(EInfo, "Block: near = (%.6f, %.6f, %.6f), far = (%.6f, %.6f, %.6f)",
					rayBlock.o.x, rayBlock.o.y, rayBlock.o.z,
					rayBlock(blockLength).x, rayBlock(blockLength).y, rayBlock(blockLength).z);
				Log(EInfo, "rayBlockDir = (%.6f, %.6f, %.6f)", 
					rayBlock.d.x, rayBlock.d.y, rayBlock.d.z);
#endif

				Float tBlock;
				PatchIntersectionRecord tempBlock;

				if (m_block->rayIntersect(rayBlock, rayBlock.mint, rayBlock.maxt,
					tBlock, (void*)(&tempBlock))) {
					// intersection in block space
					Intersection itsBlock;
					m_block->fillIntersectionRecord(rayBlock, (void*)(&tempBlock), itsBlock);

					Float tWorld = Epsilon + (maxT - minT - 2.0 * Epsilon) * (tBlock / blockLength);
					//Float tWorld = tBlock / blockLength * (maxT - minT);
					tTemp += tWorld;

					t = tTemp;
					Point itsPWorld = ray(t);
					Point4 bb;

					if (!m_shell.m_tetra[id].inside(m_shell.m_vtxPosition, itsPWorld, bb)) {
						//Log(EError, "should not happen!");
						break;
					}

					if (tmp) {
						ShellmapIntersectionRecord &temp = *((ShellmapIntersectionRecord*)tmp);
						temp.normBlock = itsBlock.geoFrame.n;
						temp.dpduBlock = itsBlock.dpdu;
						temp.dpdvBlock = itsBlock.dpdv;
						temp.bb = bb;
						temp.tetraId = id;
					}

#ifdef SHELLMAP_HEIGHTFIELD_DEBUG
					Log(EInfo, "*** found intersection ***");
					Log(EInfo, "blockT = %.8f", tBlock);
					Log(EInfo, "block.p = (%.6f, %.6f, %.6f)", itsBlock.p.x, itsBlock.p.y, itsBlock.p.z);
					Log(EInfo, "final dist = %.8f", t);
					Log(EInfo, "its.p = (%.6f, %.6f, %.6f)", itsPWorld.x, itsPWorld.y, itsPWorld.z);
#endif

					return true;
				}
			}
            
            if (maxT >= _maxt)
                break;

            // move to the next tetrahedron
            id = m_shell.m_link[4 * id + maxF];
                
            if (id < 0) {
                // ray goes out of shellmap, but could intersect with it again
                Intersection its;
                Ray rayTemp(ray, maxT + 2.0 * Epsilon, std::numeric_limits<Float>::infinity());
				if (!m_shell.m_btree->rayIntersect(rayTemp, its) || its.t > _maxt) {
#ifdef SHELLMAP_HEIGHTFIELD_DEBUG
					Log(EInfo, "Will not intersect the shellmap again.");
#endif
					break;
				}
            
                if (!m_shell.lookupPoint(ray(its.t + Epsilon), id0))
                    break;
                id = static_cast<int>(id0);
            }
                
            if (!m_shell.m_tetra[id].rayIntersect(m_shell.m_vtxPosition, _ray, minT, minF, maxT, maxF))
                break;
            tTemp = minT;
        }

        t = std::numeric_limits<Float>::infinity();
        return false;
    }

    bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
        Float t;
        return rayIntersect(ray, mint, maxt, t, NULL);
    }

    void fillIntersectionRecord(const Ray &ray, const void *_tmp, Intersection &its) const {
        ShellmapIntersectionRecord &temp = *((ShellmapIntersectionRecord*)_tmp);
        
        its.p = ray(its.t);

        int id = temp.tetraId;
        const uint32_t *tmp = m_shell.m_tetra[id].idx;
        Point tex;
        Vector norm;
        TangentSpace tang;

        tex = m_shell.m_vtxTexcoord[tmp[0]]*temp.bb.x +
            m_shell.m_vtxTexcoord[tmp[1]]*temp.bb.y +
            m_shell.m_vtxTexcoord[tmp[2]]*temp.bb.z +
            m_shell.m_vtxTexcoord[tmp[3]]*temp.bb.w;
        norm = m_shell.m_vtxNormal[tmp[0]]*temp.bb.x +
            m_shell.m_vtxNormal[tmp[1]]*temp.bb.y +
            m_shell.m_vtxNormal[tmp[2]]*temp.bb.z +
            m_shell.m_vtxNormal[tmp[3]]*temp.bb.w;
        tang.dpdu = m_shell.m_vtxTangent[tmp[0]].dpdu*temp.bb.x +
            m_shell.m_vtxTangent[tmp[1]].dpdu*temp.bb.y +
            m_shell.m_vtxTangent[tmp[2]].dpdu*temp.bb.z +
            m_shell.m_vtxTangent[tmp[3]].dpdu*temp.bb.w;
        tang.dpdv = m_shell.m_vtxTangent[tmp[0]].dpdv*temp.bb.x +
            m_shell.m_vtxTangent[tmp[1]].dpdv*temp.bb.y +
            m_shell.m_vtxTangent[tmp[2]].dpdv*temp.bb.z +
            m_shell.m_vtxTangent[tmp[3]].dpdv*temp.bb.w;

#ifdef SHELLMAP_HEIGHTFIELD_DEBUG
        Log(EInfo, "===== shellmap frame =====");
        Log(EInfo, "n = (%.6f, %.6f, %.6f)", norm.x, norm.y, norm.z);
        Log(EInfo, "s = (%.6f, %.6f, %.6f)", tang.dpdu.x, tang.dpdu.y, tang.dpdu.z);
        Log(EInfo, "t = (%.6f, %.6f, %.6f)", tang.dpdv.x, tang.dpdv.y, tang.dpdv.z);
        Log(EInfo, "===== heightmap space =====");
#endif

        its.uv = Point2(tex.x, tex.y);

        Vector dpduTex = temp.dpduBlock;
        Vector dpdvTex = temp.dpdvBlock;
        Normal normTex = temp.normBlock;

// 		Vector dpduTex = m_dataToTexture(temp.dpduBlock);
// 		Vector dpdvTex = m_dataToTexture(temp.dpdvBlock);
// 		Normal normTex = normalize(m_dataToTexture(temp.normBlock));
       
        Vector dpduWorld = dpduTex.x * tang.dpdu + dpduTex.y * tang.dpdv + dpduTex.z * norm;
        Vector dpdvWorld = dpdvTex.x * tang.dpdu + dpdvTex.y * tang.dpdv + dpdvTex.z * norm;
        Normal normWorld = normTex.x * tang.dpdu + normTex.y * tang.dpdv + normTex.z * norm;
        normWorld = normalize(normWorld);

        its.dpdu = dpduWorld;
        its.dpdv = dpdvWorld;
        
#ifdef SHELLMAP_HEIGHTFIELD_DEBUG
        Log(EInfo, "===== normal transform =====");
        Log(EInfo, "nBlock = (%.6f, %.6f, %.6f)", temp.normBlock.x, temp.normBlock.y, temp.normBlock.z);
        Log(EInfo, "nTex = (%.6f, %.6f, %.6f)", normTex.x, normTex.y, normTex.z);
        Log(EInfo, "nWorld = (%.6f, %.6f, %.6f)", normWorld.x, normWorld.y, normWorld.z);
#endif

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
        its.instance = NULL;
		its.hasUVPartials = false;
        its.time = ray.time;

		its.baseFrame.s = tang.dpdu;
		its.baseFrame.t = normalize(tang.dpdv - dot(tang.dpdv, its.baseFrame.s) * its.baseFrame.s);
		its.baseFrame.n = cross(its.baseFrame.s, its.baseFrame.t);

		// compute differential geometry
		if (m_useMacroDeform) {
			Float step = 1e-4;
			Point px, py, tx, ty, bo, bx, by;
			bo = m_textureToData.transformAffine(tex);

			px = its.p + its.baseFrame.s * step;
			if (m_shell.lookupPoint(px, tx)) {
				bx = m_textureToData.transformAffine(tx);
				its.dudx = (bx.x - bo.x) / step;
				its.dvdx = (bx.y - bo.y) / step;
			} else {
				px = its.p - its.baseFrame.s * step;
				if (m_shell.lookupPoint(px, tx)) {
					bx = m_textureToData.transformAffine(tx);
					its.dudx = (bo.x - bx.x) / step;
					its.dvdx = (bo.y - bx.y) / step;
				} else {
					its.dudx = 1.0; its.dvdx = 0.0;
				}
			}

			py = its.p + its.baseFrame.t * step;
			if (m_shell.lookupPoint(py, ty)) {
				by = m_textureToData.transformAffine(ty);
				its.dudy = (by.x - bo.x) / step;
				its.dvdy = (by.y - bo.y) / step;
			} else {
				py = its.p - its.baseFrame.t * step;
				if (m_shell.lookupPoint(py, ty)) {
					its.dudy = (bo.x - by.x) / step;
					its.dvdy = (bo.y - by.y) / step;
				} else {
					its.dudy = 0.0; its.dvdy = 1.0;
				}
			}
		}
    }

    //void getNormalDerivative(const Intersection &its, Vector &dndu, Vector &dndv, 
    //    bool shadingFrame) const {
    //}

    void configure() {
        if (!m_ready) {
            Shape::configure();

            if (m_block.get() == NULL) 
                Log(EError, "No embedded heightfield specified!");

            m_worldToObject = m_objectToWorld.inverse();
            
            AABB blockAABB = m_block->getAABB();
            Float maxz = std::max(blockAABB.max.z, m_maxHeight);
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
            Log(EInfo, "Building a BVH for the shell mesh...");
            m_shell.configure(m_objectToWorld);
            Log(EInfo, "Shell mesh loaded: %u tetrahedra, tree depth: %u",
                m_shell.getTetrahedronCount(), m_shell.getTreeDepth());

            m_aabb = m_shell.getAABB();
			m_tetrahedronCount = m_shell.getTetrahedronCount();
                
            m_ready = true;

            Log(EInfo, "%s", toString().c_str());
			Log(EInfo, "shell area = %.6f", m_shell.getSurfaceArea());
        }
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Shape)) && name == "baseHeightfield") {
            Assert(m_block == NULL);
            m_block = static_cast<Shape*>(child);
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
        if (p.z < Epsilon || p.z > 1.0 - Epsilon) {
            //Log(EError, "bad z = %.8f", p.z);
            p.z = math::clamp(p.z, 0.0, 1.0);
        }
        p.x -= math::floorToInt(p.x);
        p.y -= math::floorToInt(p.y);
    }

protected:
    std::string m_shellFilename;
    ref<Shape> m_block;
    
    bool m_ready;
    TetrahedronMesh m_shell;
    Float m_maxHeight;
    Transform m_worldToObject, m_objectToWorld;
    Transform m_textureToData, m_dataToTexture;
    
	AABB m_aabb;
	int m_tetrahedronCount;

	bool m_useMacroDeform;
};

MTS_IMPLEMENT_CLASS_S(ShellmapHeightfield, false, Shape)
MTS_EXPORT_PLUGIN(ShellmapHeightfield, "Shellmapped height field intersection shape");
MTS_NAMESPACE_END
