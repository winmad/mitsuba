#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <boost/filesystem/path.hpp>

MTS_NAMESPACE_BEGIN

class samplePatchVis : public Utility {
public:
	int run(int argc, char **argv) {
		m_sqrtSpp = std::atoi(argv[2]);
		Float umin = std::atof(argv[3]);
		Float umax = std::atof(argv[4]);
		Float vmin = std::atof(argv[5]);
		Float vmax = std::atof(argv[6]);
		//m_shadowOption = std::atoi(argv[10]);
		
		ParameterMap params;
		params["heightFilename"] = std::string(argv[7]);
		params["xscale"] = "1";
		params["yscale"] = "1";
		params["xoffset"] = "0";
		params["yoffset"] = "0";

		m_scene = loadScene(argv[1], params);

		m_aabb = AABB2(Point2(umin, vmin), Point2(umax, vmax));
		m_scene->initialize();
		//m_wi = normalize(m_wi);
		//m_wo = normalize(m_wo);
		m_spp = m_sqrtSpp * m_sqrtSpp;

		Properties props = Properties("independent");
		props.setInteger("seed", 19931004);
		m_sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), props));
		m_sampler->configure();

		// generate samples
		std::vector<std::vector<Point2> > directionSamples(m_sqrtSpp, std::vector<Point2>(m_sqrtSpp));
		std::vector<std::vector<Point2> > positionSamples(m_sqrtSpp, std::vector<Point2>(m_sqrtSpp));
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				directionSamples[i][j] = m_sampler->next2D();
				positionSamples[i][j] = m_sampler->next2D();
			}
		}

		std::vector<std::vector<Vector> > wi(m_sqrtSpp, std::vector<Vector>(m_sqrtSpp));
		std::vector<std::vector<double> > pdfWi(m_sqrtSpp, std::vector<double>(m_sqrtSpp));
		std::vector<std::vector<Vector> > wo(m_sqrtSpp, std::vector<Vector>(m_sqrtSpp));
		std::vector<std::vector<double> > pdfWo(m_sqrtSpp, std::vector<double>(m_sqrtSpp));

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
				// sample incoming direction wi
				Vector wiLocal = warp::squareToCosineHemisphere(directionSamples[i][j]);
				pdfWi[i][j] = warp::squareToCosineHemispherePdf(wiLocal);

				Ray ray;
				bool frontFaced = true;
				
				if (useHeightfield)
					frontFaced = sampleRayUV(i, j, wiLocal, positionSamples[i][j], 
						ray, normals[i][j]);
				wi[i][j] = ray.d;

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

		double avgShadowingMasking = 0.0;

		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				avgShadowingMasking += vis[i][j][2] * std::max(0.0, dot(wi[i][j], normals[i][j])) / pdfWi[i][j];
			}
		}
		avgShadowingMasking /= (double)m_spp;
		
		Log(EInfo, "===========================================");
		Log(EInfo, "Average shadowing masking: %.8f", avgShadowingMasking);

		FILE *fp = fopen("tmp_result.txt", "w");
		fprintf(fp, "%.8f", avgShadowingMasking);
		fprintf(fp, "\n");
		fclose(fp);
		
		return 0;
	}

	bool sampleRayUV(int i, int j, const Vector &wiLocal, const Point2 &sample, 
		Ray &ray, Normal &normal) {
		Point2 uv;
		uv.x = m_aabb.min.x + (j + sample.x) / (double)m_sqrtSpp * (m_aabb.max.x - m_aabb.min.x);
		uv.y = m_aabb.min.y + (i + sample.y) / (double)m_sqrtSpp * (m_aabb.max.y - m_aabb.min.y);

		//uv.x = m_aabb.min.x + (j + 0.5) / (double)m_sqrtSpp * (m_aabb.max.x - m_aabb.min.x);
		//uv.y = m_aabb.min.y + (i + 0.5) / (double)m_sqrtSpp * (m_aabb.max.y - m_aabb.min.y);

		Point pos;
		Normal norm;
		m_hmap->getPosAndNormal(uv, &pos, &norm);

		Frame frame(norm);
		Vector wi = frame.toWorld(wiLocal);
		
		normal = norm;
		ray = Ray(pos + wi * ShadowEpsilon, wi, 0);

		return true;
	}

	ref<Scene> m_scene;
	Shape *m_hmap;
	ref<Sampler> m_sampler;
	int m_sqrtSpp, m_spp;
	int m_shadowOption;
	//Vector m_wi;
	Vector m_wo;
	AABB2 m_aabb;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(samplePatchVis, "Compute the average visibility given a patch in the heightmap")
MTS_NAMESPACE_END
