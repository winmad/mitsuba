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

class CalcPatchVisibility : public Utility {
public:
	int run(int argc, char **argv) {
		m_scene = loadScene(argv[1]);
		m_size = std::atoi(argv[2]);
		m_sqrtSpp = std::atoi(argv[3]);
		Float xmin = std::atof(argv[4]);
		Float xmax = std::atof(argv[5]);
		Float ymin = std::atof(argv[6]);
		Float ymax = std::atof(argv[7]);
		m_shadowOption = std::atoi(argv[8]);
		fs::path filename(argv[9]);

		m_aabb = AABB2(Point2(xmin, ymin), Point2(xmax, ymax));
		m_scene->initialize();
		m_spp = m_sqrtSpp * m_sqrtSpp;

		Properties props = Properties("independent");
		props.setInteger("seed", 19931004);
		m_sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), props));
		m_sampler->configure();

		m_res = new Bitmap(Bitmap::ELuminance, Bitmap::EFloat32, Vector2i(m_size, m_size));
		float *data = m_res->getFloat32Data();

		m_hmap = m_scene->getShapes()[0];

		/*
		int m = 4;
		for (int i = 0; i <= m; i++) {
			for (int j = 0; j <= m; j++) {
				double x = m_aabb.min.x + (double)j / (double)m * (m_aabb.max.x - m_aabb.min.x);
				double y = m_aabb.min.y + (double)i / (double)m * (m_aabb.max.y - m_aabb.min.y);
				Point o(x, y, 1e2);
				Ray ray(o, Vector(0, 0, -1), 0);
				Intersection its;
				m_scene->rayIntersect(ray, its);
				Float h = hmap->getHeight(o);
				if (std::abs(h - its.p.z) > Epsilon) {
					Log(EInfo, "different heights: %.6f %.6f", h, its.p.z);
				}

				Vector normal = hmap->getNormal(o);
				Vector diff = normal - its.geoFrame.n;
				if (diff.length() > Epsilon) {
					Log(EInfo, "different normals: (%.6f, %.6f, %.6f), (%.6f, %.6f, %.6f)",
						normal.x, normal.y, normal.z,
						its.geoFrame.n.x, its.geoFrame.n.y, its.geoFrame.n.z);
				}
			}
		}
		return 0;
		*/

#pragma omp parallel for
		for (int r = 0; r < m_size; r++) {
			for (int c = 0; c < m_size; c++) {
				double res = 0.0f;
				for (int i = 0; i < m_sqrtSpp; i++) {
					for (int j = 0; j < m_sqrtSpp; j++) {
						double y = (r + 0.5) / (double)m_size * 2.0 - 1.0;
						double x = (c + 0.5) / (double)m_size * 2.0 - 1.0;
						double z = 1.0 - x * x - y * y;
						if (z < Epsilon)
							continue;
						z = std::sqrt(z);
						Vector wo(x, y, z);

						Point o = sampleRayOrigin(i, j, m_sampler);
						Ray ray(o + wo * Epsilon, wo, 0);
						//Ray ray(o, wo, 0);
						Ray invRay(o + wo * 1e3, -wo, 0);
						//Log(EInfo, "%.6f, %.6f, %.6f", wo.x, wo.y, wo.z);

						Intersection its;
						double vis = 0.0f;

						// verify ray.o is above the heightfield
						bool isInside = false;
						Float h = m_hmap->getHeight(ray.o);
						if (h > ray.o.z) {
							isInside = true;
							//Log(EInfo, "ray origin is inside, %.6f, %.6f", ray.o.z, its.p.z);
						}

						if (m_shadowOption == 1) {
							bool flag = m_scene->rayIntersect(ray, its);
							if (!isInside && (!flag || (its.isValid() && !m_aabb.contains(Point2(its.p.x, its.p.y))))) {
								vis = 1.0f;
// 								if (its.isValid() && dot(its.geoFrame.n, -ray.d) < Epsilon) {
// 									vis = 0.f;
// 									Log(EInfo, "back face normal = (%.6f, %.6f, %.6f)", its.geoFrame.n.x,
// 										its.geoFrame.n.y, its.geoFrame.n.z);
// 								}
							}
						}
						else if (m_shadowOption == 2) {
							m_scene->rayIntersect(invRay, its);
							if (its.isValid() && m_aabb.contains(Point2(its.p.x, its.p.y)))
								vis = 1.0f;
						}
						else if (m_shadowOption == 3) {
							//m_scene->rayIntersect(invRay, its);
							//if (!isInside && its.isValid() && (its.p - o).length() < Epsilon) {
							//	vis = 1.0f;
							//}

							bool flag = m_scene->rayIntersect(ray);
							if (!isInside && !flag) {
								vis = 1.0f;
							}
						}

						res += vis;
					}
				}

				res /= (double)m_spp;
				data[(m_size - r - 1) * m_size + c] = res;
			}
		}

		Log(EInfo, "Finish computation.");

		double avgVis = 0.0f;
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

		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		m_res->write(Bitmap::EOpenEXR, stream);

		return 0;
	}

	Point sampleRayOrigin(int i, int j, Sampler *sampler) {
		double x = m_aabb.min.x + (j + sampler->next1D()) / (double)m_sqrtSpp * (m_aabb.max.x - m_aabb.min.x);
		double y = m_aabb.min.y + (i + sampler->next1D()) / (double)m_sqrtSpp * (m_aabb.max.y - m_aabb.min.y);
		Point o(x, y, 0);
		double z = m_hmap->getHeight(o);
		return Point(x, y, z);

// 		Point o(x, y, 1e2);
// 		Ray ray(o, Vector(0, 0, -1.0f), 0);
// 		Intersection its;
// 		m_scene->rayIntersect(ray, its);
//		return its.p;
	}

	ref<Scene> m_scene;
	Shape *m_hmap;
	ref<Sampler> m_sampler;
	int m_sqrtSpp, m_spp;
	int m_size;
	AABB2 m_aabb;
	int m_shadowOption;
	ref<Bitmap> m_res;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(CalcPatchVisibility, "Compute the spherical visibility function of a patch in the heightmap")
MTS_NAMESPACE_END
