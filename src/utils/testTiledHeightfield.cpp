#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <boost/filesystem/path.hpp>

MTS_NAMESPACE_BEGIN

namespace {
    struct PatchIntersectionRecord {
        Point p;
        int x, y;
        int blockX, blockY;
    };
}

class TestTiledHeightfield : public Utility {
public:
	int run(int argc, char **argv) {
        m_scene = loadScene(argv[1]);
        Point o;
        o.x = std::atof(argv[2]);
        o.y = std::atof(argv[3]);
        o.z = std::atof(argv[4]);
        Vector d;
        d.x = std::atof(argv[5]);
        d.y = std::atof(argv[6]);
        d.z = std::atof(argv[7]);
        d = normalize(d);
		
		Properties props = Properties("independent");
		props.setInteger("seed", 19931004);
		m_sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), props));
		m_sampler->configure();

        m_tiledHmap = m_scene->getShapes()[0];
        m_hmap = m_scene->getShapes()[1];

        compare(o, d, true);

        for (int i = 0; i < 0; i++) {
            o.x = (m_sampler->next1D() * 2.0 - 1) * 500.0;
            o.y = (m_sampler->next1D() * 2.0 - 1) * 500.0;
            o.z = 20;

            d.x = m_sampler->next1D();
            d.y = m_sampler->next1D();
            d.z = -1;
            d = normalize(d);

            if (!compare(o, d)) {
                Log(EInfo, "--- fail cmp test %d ---", i);
                compare(o, d, true);
            }
        }

        bool flag = true;
        for (int i = 0; i < 5; i++) {
            o.x = (m_sampler->next1D() * 2.0 - 1) * 1024.0;
            o.y = (m_sampler->next1D() * 2.0 - 1) * 1024.0;
            o.z = 20;

            if (!testGetHeightAndNormal(o, true)) {
                flag = false;
                Log(EInfo, "--- fail unit test %d ---", i);
                testGetHeightAndNormal(o, true);
            }
        }
        if (flag)
            Log(EInfo, "--- unit test success ---");
		return 0;
    }
    
    bool compare(const Point &o, const Vector &d, bool verbose = false) {
        Ray ray(o, d, 0);
        Float mint = 0;
        Float maxt = std::numeric_limits<Float>::infinity();
        
        Float t1;
        PatchIntersectionRecord temp1;
        m_tiledHmap->rayIntersect(ray, mint, maxt, t1, (void*)(&temp1));
        if (verbose) {
            Log(EInfo, "===== TiledHeightfield =====");
            Log(EInfo, "t = %.6f", t1);
            Log(EInfo, "point = (%.6f, %.6f, %.6f)", temp1.p.x, temp1.p.y, temp1.p.z);
            Log(EInfo, "x = %d, y = %d", temp1.x, temp1.y);
            Log(EInfo, "blockX = %d, blockY = %d", temp1.blockX, temp1.blockY);
        }
        Intersection its1;
        m_tiledHmap->fillIntersectionRecord(ray, (void*)(&temp1), its1);
        if (verbose) {
            Log(EInfo, "p = (%.6f, %.6f, %.6f)", its1.p.x, its1.p.y, its1.p.z);
            Log(EInfo, "n = (%.6f, %.6f, %.6f)", its1.geoFrame.n.x, its1.geoFrame.n.y, its1.geoFrame.n.z);
        }

        Float t2;
        PatchIntersectionRecord temp2;
        m_hmap->rayIntersect(ray, mint, maxt, t2, (void*)(&temp2));
        if (verbose) {
            Log(EInfo, "===== Heightfield =====");
            Log(EInfo, "t = %.6f", t2);
            Log(EInfo, "point = (%.6f, %.6f, %.6f)", temp2.p.x, temp2.p.y, temp2.p.z);
            Log(EInfo, "x = %d, y = %d", temp2.x, temp2.y);
        }
        Intersection its2;
        m_hmap->fillIntersectionRecord(ray, (void*)(&temp2), its2);
        if (verbose) {
            Log(EInfo, "p = (%.6f, %.6f, %.6f)", its2.p.x, its2.p.y, its2.p.z);
            Log(EInfo, "n = (%.6f, %.6f, %.6f)", its2.geoFrame.n.x, its2.geoFrame.n.y, its2.geoFrame.n.z);
        }

        bool flag = true;
        if ((its1.p - its2.p).length() > 1e-3f)
            flag = false;
        if ((its1.geoFrame.n - its2.geoFrame.n).length() > 1e-3f)
            flag = false;
        return flag;
    }

    bool testGetHeightAndNormal(const Point &o, bool verbose = false) {
        Ray ray(o, Vector(0, 0, -1), 0);
        Float mint = 0;
        Float maxt = std::numeric_limits<Float>::infinity();
        Float t;
        PatchIntersectionRecord temp;
        m_tiledHmap->rayIntersect(ray, mint, maxt, t, (void*)(&temp));
        Intersection its;
        m_tiledHmap->fillIntersectionRecord(ray, (void*)(&temp), its);

        Float h = m_tiledHmap->getHeight(o);
        Vector n = m_tiledHmap->getNormal(o);

        if (verbose) {
            Log(EInfo, "h: %.6f, %.6f", h, its.p.z);
            Log(EInfo, "n: (%.6f, %.6f, %.6f), (%.6f, %.6f, %.6f)", n.x, n.y, n.z,
                its.geoFrame.n.x, its.geoFrame.n.y, its.geoFrame.n.z);
        }

        bool flag = true;
        if (std::abs(h - its.p.z) > 1e-3f)
            flag = false;
        if ((n - its.geoFrame.n).length() > 1e-3f)
            flag = false;
        return flag;
    }

	ref<Scene> m_scene;
    Shape *m_hmap;
    Shape *m_tiledHmap;
	ref<Sampler> m_sampler;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(TestTiledHeightfield, "Test tiled heightfield")
MTS_NAMESPACE_END
