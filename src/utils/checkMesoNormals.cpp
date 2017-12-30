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
#include <mitsuba/render/texture.h>
#include <boost/filesystem/path.hpp>

MTS_NAMESPACE_BEGIN

class CheckMesoNormals : public Utility {
public:
	typedef std::vector<std::vector<Normal> > Normal2DArray;

	int run(int argc, char **argv) {
		// avoid negative value
		m_offset = 5.0f;

		ParameterMap params;
		params["heightmap"] = std::string(argv[1]);
		m_scene = loadScene(fs::path("heightmap_LEADR.xml"), params);

		Properties props;
		m_scale = std::atoi(argv[3]);
		m_sqrtSpp = std::atoi(argv[4]);
		m_size = std::atoi(argv[5]);

		char filename[256];
		sprintf(filename, "%s/moments0_%dx.exr", argv[2], m_scale);
		props = Properties("bitmap");
		props.setString("filename", filename);
		props.setFloat("gamma", 1.0);
		m_moments0 = static_cast<Texture *>(PluginManager::getInstance()->
			createObject(MTS_CLASS(Texture), props));
		m_moments0->configure();

		sprintf(filename, "%s/moments1_%dx.exr", argv[2], m_scale);
		props = Properties("bitmap");
		props.setString("filename", filename);
		props.setFloat("gamma", 1.0);
		m_moments1 = static_cast<Texture *>(PluginManager::getInstance()->
			createObject(MTS_CLASS(Texture), props));
		m_moments1->configure();

		m_scene->initialize();
		m_spp = m_sqrtSpp * m_sqrtSpp;

		props = Properties("independent");
		props.setInteger("seed", 19931004);
		m_sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), props));
		m_sampler->configure();

		// generate samples
		std::vector<std::vector<Point2> > positionSamples(m_sqrtSpp, std::vector<Point2>(m_sqrtSpp));
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				positionSamples[i][j] = m_sampler->next2D();
			}
		}

		m_hmap = m_scene->getShapes()[0];
		bool useHeightfield = true;
		if (m_hmap->getClass()->getName() != "Heightfield" &&
			m_hmap->getClass()->getName() != "TiledHeightfield" &&
			m_hmap->getClass()->getName() != "ShellmapHeightfield") {
			Log(EInfo, "%s", m_hmap->getClass()->getName().c_str());
			useHeightfield = false;
		}

		std::vector<std::vector<Normal> > geoNormals(m_size, std::vector<Normal>(m_size));
		std::vector<std::vector<Normal> > leadrNormals(m_size, std::vector<Normal>(m_size));

#pragma omp parallel for
		for (int y = 0; y < m_size; y++) {
			for (int x = 0; x < m_size; x++) {
				Spectrum moments0(0.0f);
				Spectrum moments1(0.0f);

				geoNormals[y][x] = Normal(0.0f);
				leadrNormals[y][x] = Normal(0.0f);

				for (int i = 0; i < m_sqrtSpp; i++) {
					for (int j = 0; j < m_sqrtSpp; j++) {
						Point2 uv;
						uv.x = (x + (j + positionSamples[i][j].x) / (double)m_sqrtSpp) / (double)m_size;
						uv.y = (y + (i + positionSamples[i][j].y) / (double)m_sqrtSpp) / (double)m_size;

						Point pos;
						Normal norm;
						m_hmap->getPosAndNormal(uv, &pos, &norm);
						//Log(EInfo, "(%.6f, %.6f, %.6f)", norm.x, norm.y, norm.z);
						geoNormals[y][x] += norm;

						Point o(pos.x, pos.y, 1e2);
						Ray ray(o, Vector(0, 0, -1.0f), 0);
						Intersection its;
						m_scene->rayIntersect(ray, its);

						moments0 += m_moments0->eval(its) - Spectrum(m_offset);
						moments1 += m_moments1->eval(its) - Spectrum(m_offset);
						
					}
				}

				geoNormals[y][x] = normalize(geoNormals[y][x]);

				moments0 /= (double)m_spp;
				moments1 /= (double)m_spp;
				leadrNormals[y][x] = Normal(-moments0[0], -moments0[1], 1);
				leadrNormals[y][x] = normalize(leadrNormals[y][x]);
			}
		}

		outputNormals(geoNormals, leadrNormals);
		return 0;
	}

	void outputNormals(const Normal2DArray &geoNormals, const Normal2DArray &leadrNormals) {
		int resSize = m_size;
		ref<Bitmap> res0 = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, Vector2i(resSize, resSize));
		ref<Bitmap> res1 = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, Vector2i(resSize, resSize));
		res0->clear();
		res1->clear();
		float *data0 = res0->getFloat32Data();
		float *data1 = res1->getFloat32Data();

		for (int i = 0; i < resSize; i++) {
			for (int j = 0; j < resSize; j++) {
				int idx = resSize * i + j;
				for (int c = 0; c < 3; c++) {
					data0[3 * idx + c] = (geoNormals[i][j][c] + 1.0) * 0.5;
					data1[3 * idx + c] = (leadrNormals[i][j][c] + 1.0) * 0.5;
				}
			}
		}

		Log(EInfo, "Output normals...");

		char filename[256];
		sprintf(filename, "geoNormals_%dx.exr", m_scale);
		ref<FileStream> stream = new FileStream(fs::path(filename), FileStream::ETruncWrite);
		res0->write(Bitmap::EOpenEXR, stream);

		sprintf(filename, "leadrNormals_%dx.exr", m_scale);
		stream = new FileStream(fs::path(filename), FileStream::ETruncWrite);
		res1->write(Bitmap::EOpenEXR, stream);

		Log(EInfo, "Output finish");
	}

	ref<Scene> m_scene;
	Shape *m_hmap;
	ref<Texture> m_moments0;
	ref<Texture> m_moments1;
	ref<Sampler> m_sampler;
	int m_sqrtSpp, m_spp;
	int m_size;
	int m_scale;
	Float m_offset;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(CheckMesoNormals, "Output meso-normals")
MTS_NAMESPACE_END
