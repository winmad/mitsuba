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

class CalcPatchDirVisSimple : public Utility {
public:
	int run(int argc, char **argv) {
		m_wi = Vector(std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]));
		m_sqrtSpp = std::atoi(argv[5]);
		Float xmin = std::atof(argv[6]);
		Float xmax = std::atof(argv[7]);
		Float ymin = std::atof(argv[8]);
		Float ymax = std::atof(argv[9]);
		ParameterMap params;
		if (argc > 10) {
			params["height"] = std::string(argv[10]);
			params["xyscale"] = "16";
		}
		m_scene = loadScene(argv[1], params);

		m_aabb = AABB2(Point2(xmin, ymin), Point2(xmax, ymax));
		m_scene->initialize();
		m_wi = normalize(m_wi);
		m_spp = m_sqrtSpp * m_sqrtSpp;

		Properties props = Properties("independent");
		props.setInteger("seed", 19931004);
		m_sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), props));
		m_sampler->configure();

		std::vector<std::vector<double> > vis(m_sqrtSpp, std::vector<double>(m_sqrtSpp));
		std::vector<std::vector<Normal> > normals(m_sqrtSpp, std::vector<Normal>(m_sqrtSpp));
		std::vector<std::vector<double> > weights(m_sqrtSpp, std::vector<double>(m_sqrtSpp));
		m_hmap = m_scene->getShapes()[0];
		bool useHeightfield = true;
		if (m_hmap->getClass()->getName() != "Heightfield") {
			Log(EInfo, "%s", m_hmap->getClass()->getName().c_str());
			useHeightfield = false;
		}

#pragma omp parallel for
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				Point o = sampleRayOrigin(i, j, m_sampler, useHeightfield);
				Ray ray(o + m_wi * Epsilon, m_wi, 0);
				double tmp = 0.0f;

				// verify ray.o is above the heightfield
				bool isInside = false;

				if (useHeightfield) {
					Float h = m_hmap->getHeight(ray.o);
					if (h > ray.o.z) {
						isInside = true;
						//Log(EInfo, "ray origin is inside, %.6f, %.6f", ray.o.z, its.p.z);
					}
				}
				else {
					Intersection its;
					Ray vRay(Point(ray.o.x, ray.o.y, 1e2), Vector(0, 0, -1), 0);
					m_scene->rayIntersect(vRay, its);
					if (its.p.z > ray.o.z) {
						isInside = true;
					}
				}

				Normal normal;
				//normal = m_hmap->getNormal(o);
				Ray vRay(Point(o.x, o.y, 1e2), Vector(0, 0, -1), 0);
				Intersection vIts;
				m_scene->rayIntersect(vRay, vIts);
				normal = vIts.shFrame.n;
				//normal = vIts.geoFrame.n;
				normals[i][j] = normal;
				if (normal.z > 0) {
					weights[i][j] = 1.0 / normal.z / (double)m_spp;
				}
				else {
					weights[i][j] = 0.0;
				}

				// local visibility
				Intersection its;
				bool flag = m_scene->rayIntersect(ray, its);
				if (!isInside && dot(m_wi, normal) > 0 &&
					(!flag || (its.isValid() && !m_aabb.contains(Point2(its.p.x, its.p.y))))) {
					tmp = 1.0; //std::max(0.f, dot(m_wi, normal));
				}

				vis[i][j] = tmp;
			}
		}
		Log(EInfo, "Finish sampling.");

		double avgVis = 0.0;
		double totWeight = 0.0;

		double avgVisCos = 0.0;
		Normal avgNormal(0.0);
		double avgProjArea = 0.0;

		double totVisPosProjArea = 0.0;
		double totPosProjArea = 0.0;

		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				double cosM = dot(m_wi, normals[i][j]);
				if (cosM > 0) {
					avgVis += vis[i][j] * cosM;
					totWeight += cosM;
				}

				avgVisCos += vis[i][j] * std::max(0.0, dot(m_wi, normals[i][j])) / (double)m_spp;
				avgNormal += normals[i][j] * weights[i][j];
// 				if (normals[i][j].z < 0) {
// 					Log(EInfo, "n = (%.6f, %.6f, %.6f)", normals[i][j].x, normals[i][j].y, normals[i][j].z);
// 				}
				avgProjArea += std::max(0.0, normals[i][j].z) * weights[i][j];

				totVisPosProjArea += vis[i][j] * std::max(0.0, dot(m_wi, normals[i][j])) * weights[i][j];
				totPosProjArea += std::max(0.0, dot(m_wi, normals[i][j])) * weights[i][j];
			}
		}
		avgVis /= totWeight;
		avgNormal = normalize(avgNormal);
		
		Log(EInfo, "Average normal (%.6f, %.6f, %.6f), cos_wo = %.6f, cos_wn = %.6f", 
			avgNormal.x, avgNormal.y, avgNormal.z, 
			dot(avgNormal, m_wi), avgProjArea);
		Log(EInfo, "Average cosine visibility: %.8f", avgVisCos);
		Log(EInfo, "Average visibility: %.8f", avgVis);
		Log(EInfo, "G1 = visArea / totPosArea: %.8f / %.8f = %.8f", totVisPosProjArea, totPosProjArea,
			totVisPosProjArea / totPosProjArea);

		FILE *fp = fopen("tmp_result.txt", "w");
		fprintf(fp, "%.8f\n", dot(avgNormal, m_wi) / totPosProjArea);
		fclose(fp);

		return 0;
	}

	Point sampleRayOrigin(int i, int j, Sampler *sampler, bool useHeightfield) {
		double x = m_aabb.min.x + (j + sampler->next1D()) / (double)m_sqrtSpp * (m_aabb.max.x - m_aabb.min.x);
		double y = m_aabb.min.y + (i + sampler->next1D()) / (double)m_sqrtSpp * (m_aabb.max.y - m_aabb.min.y);
		double z;
		Point o(x, y, 0);
		if (useHeightfield) {
			z = m_hmap->getHeight(o);
		}
		else {
			o.z = 1e2;
			Ray ray(o, Vector(0, 0, -1.0f), 0);
			Intersection its;
			m_scene->rayIntersect(ray, its);
			z = its.p.z;
		}
		return Point(x, y, z);
	}

	ref<Scene> m_scene;
	Shape *m_hmap;
	ref<Sampler> m_sampler;
	int m_sqrtSpp, m_spp;
	Vector m_wi;
	AABB2 m_aabb;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(CalcPatchDirVisSimple, "Compute the directional visibility given a patch in the heightmap")
MTS_NAMESPACE_END
