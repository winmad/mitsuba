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

class CalcPatchDirVis : public Utility {
public:
	int run(int argc, char **argv) {
		m_wi = Vector(std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]));
		m_sqrtSpp = std::atoi(argv[5]);
		Float umin = std::atof(argv[6]);
		Float umax = std::atof(argv[7]);
		Float vmin = std::atof(argv[8]);
		Float vmax = std::atof(argv[9]);
		//m_shadowOption = std::atoi(argv[10]);
		
		ParameterMap params;
		params["heightFilename"] = std::string(argv[10]);
		params["wix"] = argv[2];
		params["wiy"] = argv[3];
		params["wiz"] = argv[4];
		params["xscale"] = "1";
		params["yscale"] = "1";
		params["xoffset"] = "0";
		params["yoffset"] = "0";

		m_scene = loadScene(argv[1], params);

		m_aabb = AABB2(Point2(umin, vmin), Point2(umax, vmax));
		m_scene->initialize();
		m_wi = normalize(m_wi);
		m_spp = m_sqrtSpp * m_sqrtSpp;

		Properties props = Properties("independent");
		props.setInteger("seed", 19931004);
		m_sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), props));
		m_sampler->configure();

		//std::vector<std::vector<double> > vis(m_sqrtSpp, std::vector<double>(m_sqrtSpp));
		std::vector<std::vector<std::vector<double> > > vis(
			m_sqrtSpp, std::vector<std::vector<double> >(m_sqrtSpp, std::vector<double>(3)));
		std::vector<std::vector<Normal> > normals(m_sqrtSpp, std::vector<Normal>(m_sqrtSpp));
		std::vector<std::vector<double> > weights(m_sqrtSpp, std::vector<double>(m_sqrtSpp));
		
		m_hmap = m_scene->getShapes()[0];
		bool useHeightfield = true;
		if (m_hmap->getClass()->getName() != "Heightfield" &&
			m_hmap->getClass()->getName() != "TiledHeightfield" &&
			m_hmap->getClass()->getName() != "ShellmapHeightfield") {
			Log(EInfo, "%s", m_hmap->getClass()->getName().c_str());
			useHeightfield = false;
		}

#pragma omp parallel for
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				Ray ray;
				bool frontFaced = true;
				
				if (useHeightfield)
					frontFaced = sampleRayUV(i, j, m_sampler, ray, normals[i][j]);

				if (normals[i][j].z > 0) {
					weights[i][j] = 1.0 / normals[i][j].z;
				}
				else {
					weights[i][j] = 0.0;
				}

				if (!frontFaced) {
					for (int k = 0; k < 3; k++) {
						vis[i][j][k] = 0.0;
					}
					continue;
				}

				Intersection its;
				bool flag = m_scene->rayIntersect(ray, its);
				if (!flag) {
					for (int k = 0; k < 3; k++) {
						vis[i][j][k] = 1.0;
					}
					continue;
				}

				vis[i][j][2] = 0.0;
				Assert(its.isValid());
				if (!m_aabb.contains(its.uv)) {
					vis[i][j][0] = 1.0;
					vis[i][j][1] = 0.0;
					continue;
				}

				vis[i][j][0] = 0.0;
				bool distantVisible = true;
				while (flag) {
					if (!m_aabb.contains(its.uv)) {
						distantVisible = false;
						break;
					}
					ray.o = its.p + ray.d * ShadowEpsilon;
					flag = m_scene->rayIntersect(ray, its);
				}
				vis[i][j][1] = (distantVisible ? 1.0 : 0.0);
			}
		}
		Log(EInfo, "Finish sampling.");

		Vector avgNormal(0.0);
		double avgProjArea = 0.0;

		std::vector<double> totVisPosProjArea(3, 0.0);
		double totPosProjArea = 0.0;
		std::vector<double> G(3, 0.0);

		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				avgNormal += normals[i][j] * weights[i][j];
				avgProjArea += std::max(0.0, normals[i][j].z) * weights[i][j];
				
				totPosProjArea += std::max(0.0, dot(m_wi, normals[i][j])) * weights[i][j];
				for (int k = 0; k < 3; k++) {
					totVisPosProjArea[k] += vis[i][j][k] * std::max(0.0, dot(m_wi, normals[i][j])) * weights[i][j];
				}
			}
		}
		avgNormal = normalize(avgNormal);
		avgProjArea /= (double)m_spp;
		totPosProjArea /= (double)m_spp;
		for (int k = 0; k < 3; k++) {
			totVisPosProjArea[k] /= (double)m_spp;
			G[k] = totVisPosProjArea[k] / totPosProjArea;
		}
		
		Log(EInfo, "Average normal (%.6f, %.6f, %.6f), cos_wo = %.6f, cos_wn = %.6f", 
			avgNormal.x, avgNormal.y, avgNormal.z, 
			dot(avgNormal, m_wi), avgProjArea);
		Log(EInfo, "G_local = visArea / totPosArea: %.8f / %.8f = %.8f", 
			totVisPosProjArea[0], totPosProjArea, G[0]);
		Log(EInfo, "G_distant = visArea / totPosArea: %.8f / %.8f = %.8f", 
			totVisPosProjArea[1], totPosProjArea, G[1]);
		Log(EInfo, "G_total = visArea / totPosArea: %.8f / %.8f = %.8f", 
			totVisPosProjArea[2], totPosProjArea, G[2]);
		Log(EInfo, "Correlation: %.8f", G[0] * G[1] - G[2]);

		FILE *fp = fopen("tmp_result.txt", "w");
		for (int k = 0; k < 3; k++) {
			fprintf(fp, "%.8f ", G[k]);
		}
		fprintf(fp, "\n");
		fclose(fp);
		
		return 0;
	}

	bool sampleRayUV(int i, int j, Sampler *sampler, Ray &ray, Normal &normal) {
		Point2 uv;
		//uv.x = m_aabb.min.x + (j + sampler->next1D()) / (double)m_sqrtSpp * (m_aabb.max.x - m_aabb.min.x);
		//uv.y = m_aabb.min.y + (i + sampler->next1D()) / (double)m_sqrtSpp * (m_aabb.max.y - m_aabb.min.y);

		uv.x = m_aabb.min.x + (j + 0.5) / (double)m_sqrtSpp * (m_aabb.max.x - m_aabb.min.x);
		uv.y = m_aabb.min.y + (i + 0.5) / (double)m_sqrtSpp * (m_aabb.max.y - m_aabb.min.y);

		Point pos;
		Normal norm;
		m_hmap->getPosAndNormal(uv, &pos, &norm);
		
		normal = norm;
		ray = Ray(pos + m_wi * ShadowEpsilon, m_wi, 0);

		if (dot(m_wi, normal) < Epsilon)
			return false;
		else
			return true;
	}

	ref<Scene> m_scene;
	Shape *m_hmap;
	ref<Sampler> m_sampler;
	int m_sqrtSpp, m_spp;
	int m_shadowOption;
	Vector m_wi;
	AABB2 m_aabb;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(CalcPatchDirVis, "Compute the directional visibility given a patch in the heightmap")
MTS_NAMESPACE_END
