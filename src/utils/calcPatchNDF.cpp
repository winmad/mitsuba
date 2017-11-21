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

class CalcPatchNDF : public Utility {
public:
	int run(int argc, char **argv) {
		m_scene = loadScene(argv[1]);
		m_size = std::atoi(argv[2]);
		m_sqrtSpp = std::atoi(argv[3]);
		Float xmin = std::atof(argv[4]);
		Float xmax = std::atof(argv[5]);
		Float ymin = std::atof(argv[6]);
		Float ymax = std::atof(argv[7]);
		m_normal_stLim = std::atof(argv[8]);
		fs::path filename(argv[9]);

		m_aabb = AABB2(Point2(xmin, ymin), Point2(xmax, ymax));
		m_scene->initialize();
		m_spp = m_sqrtSpp * m_sqrtSpp;

		Properties props = Properties("independent");
		props.setInteger("seed", 19931004);
		m_sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), props));
		m_sampler->configure();

		m_res = new Bitmap(Bitmap::EPixelFormat::ELuminance, Bitmap::EFloat32, Vector2i(m_size, m_size));
		m_res->clear();
		float *data = m_res->getFloat32Data();

		m_hmap = m_scene->getShapes()[0];

		double normFactor = (double)m_size * m_size / (double)m_spp;
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				Vector normal = sampleNormal(i, j, m_sampler);
				normal /= m_normal_stLim;

				int c = (normal.x > 0.9999f ? m_size - 1 : math::floorToInt((normal.x + 1.0) * 0.5 * m_size));
				int r = (normal.y > 0.9999f ? m_size - 1 : math::floorToInt((normal.y + 1.0) * 0.5 * m_size));

				data[(m_size - r - 1) * m_size + c] += normFactor;
			}
		}
		Log(EInfo, "Finish binning.");

		/*
		double totWeight = 0.0f;
		for (int r = 0; r < m_size; r++) {
			double y = (r + 0.5) / (double)m_size * 2.0 - 1.0;
			for (int c = 0; c < m_size; c++) {
				double x = (c + 0.5) / (double)m_size * 2.0 - 1.0;
				double sinTheta2 = x * x + y * y;
				if (sinTheta2 >= 1.0)
					continue;
				double w = 1.0 / std::sqrt(1.0 - sinTheta2);
				avgVis += data[(m_size - r - 1) * m_size + c] * w;
				totWeight += w;
			}
		}
		avgVis /= totWeight;
		Log(EInfo, "Average vis = %.8f", avgVis);
		*/

		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		m_res->write(Bitmap::EOpenEXR, stream);

		return 0;
	}

	Vector sampleNormal(int i, int j, Sampler *sampler) {
		double x = m_aabb.min.x + (j + sampler->next1D()) / (double)m_sqrtSpp * (m_aabb.max.x - m_aabb.min.x);
		double y = m_aabb.min.y + (i + sampler->next1D()) / (double)m_sqrtSpp * (m_aabb.max.y - m_aabb.min.y);
// 		Point o(x, y, 0);
// 		return m_hmap->getNormal(o);

		Point o(x, y, 1e2);
		Ray ray(o, Vector(0, 0, -1.0f), 0);
		Intersection its;
		m_scene->rayIntersect(ray, its);
		return its.shFrame.n;
	}

	ref<Scene> m_scene;
	Shape *m_hmap;
	ref<Sampler> m_sampler;
	int m_sqrtSpp, m_spp;
	int m_size;
	AABB2 m_aabb;
	double m_normal_stLim;
	ref<Bitmap> m_res;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(CalcPatchNDF, "Compute the NDF given a patch in the heightmap")
MTS_NAMESPACE_END
