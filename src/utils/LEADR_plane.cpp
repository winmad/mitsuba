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

class LEADRplane : public Utility {
public:
	typedef std::vector<std::vector<Vector6> > Moment2DArray;
	typedef std::vector<std::vector<Vector3> > Vector2DArray;

	int run(int argc, char **argv) {
		// avoid negative value
		m_offset = 1e4;

		m_size = std::atoi(argv[2]);
		m_sqrtSpp = std::atoi(argv[3]);
		m_maxScale = std::atoi(argv[4]);

		ParameterMap params;
		//params["heightFilename"] = std::string(argv[7]);

		m_scene = loadScene(argv[1], params);

		m_scene->initialize();
		m_spp = m_sqrtSpp * m_sqrtSpp;

		Properties props = Properties("independent");
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

		std::vector<std::vector<Vector6> > moments(m_size, std::vector<Vector6>(m_size));
		Vector2DArray uvs(m_size, std::vector<Vector3>(m_size));

		m_hmap = m_scene->getShapes()[0];
		bool useHeightfield = true;
		if (m_hmap->getClass()->getName() != "Heightfield" &&
			m_hmap->getClass()->getName() != "TiledHeightfield" &&
			m_hmap->getClass()->getName() != "ShellmapHeightfield") {
			Log(EInfo, "%s", m_hmap->getClass()->getName().c_str());
			useHeightfield = false;
		}

#pragma omp parallel for
		for (int y = 0; y < m_size; y++) {
			for (int x = 0; x < m_size; x++) {
				moments[y][x] = Vector6(0.0);
				uvs[y][x] = Vector3(0.0);

				for (int i = 0; i < m_sqrtSpp; i++) {
					for (int j = 0; j < m_sqrtSpp; j++) {
						Point2 uv;
						uv.x = (x + (j + positionSamples[i][j].x) / (double)m_sqrtSpp) / (double)m_size;
						uv.y = (y + (i + positionSamples[i][j].y) / (double)m_sqrtSpp) / (double)m_size;

						Point pos;
						Normal norm;
						m_hmap->getPosAndNormal(uv, &pos, &norm);

						double slopeX = -norm.x / norm.z;
						double slopeY = -norm.y / norm.z;

						moments[y][x][0] += slopeX;
						moments[y][x][1] += slopeY;
						moments[y][x][2] += slopeX * slopeX;
						moments[y][x][3] += slopeY * slopeY;
						moments[y][x][4] += slopeX * slopeY;

						uvs[y][x] += Vector3(uv.x, uv.y, 0);
					}
				}

				moments[y][x] /= (double)m_spp;
				uvs[y][x] /= (double)m_spp;

				/*
				Float sigmaX2 = moments[y][x][2] - moments[y][x][0] * moments[y][x][0];
				Float sigmaY2 = moments[y][x][3] - moments[y][x][1] * moments[y][x][1];
				Float cxy = moments[y][x][4] - moments[y][x][0] * moments[y][x][1];
				if (sigmaX2 * sigmaY2 - cxy * cxy < 0) {
					Log(EInfo, "%.8f, %.8f, %.8f, %.8f, %.8f", moments[y][x][0], moments[y][x][1],
						moments[y][x][2], moments[y][x][3], moments[y][x][4]);
					//Log(EInfo, "cov matrix = (%.8f, %.8f; %.8f, %.8f)", sigmaX2, cxy, cxy, sigmaY2);
				}
				*/
			}
		}
		
		outputBitmap(moments, 1);
		outputBitmap(uvs, "uv.exr");
		momentsMipmap(moments);

		return 0;
	}

	void momentsMipmap(const Moment2DArray &moments) {
		int scale = 2;

		Moment2DArray momentsNow, momentsNext;
		momentsNow = moments;

		while (scale <= m_maxScale) {
			int sizeNext = m_size / scale;
			momentsNext = std::vector<std::vector<Vector6> >(sizeNext, std::vector<Vector6>(sizeNext));

#pragma omp parallel for
			for (int i = 0; i < sizeNext; i++) {
				for (int j = 0; j < sizeNext; j++) {
					momentsNext[i][j] = (momentsNow[2 * i][2 * j] + momentsNow[2 * i][2 * j + 1] +
						momentsNow[2 * i + 1][2 * j] + momentsNow[2 * i + 1][2 * j + 1]) * 0.25;
				}
			}

			outputBitmap(momentsNext, scale);

			scale *= 2;
			momentsNow = momentsNext;
		}
	}

	void outputBitmap(const Moment2DArray &moments, int scale) {
		int resSize = m_size / scale;
		ref<Bitmap> res0 = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, Vector2i(resSize, resSize));
		ref<Bitmap> res1 = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, Vector2i(resSize, resSize));
		res0->clear();
		res1->clear();
		float *data0 = res0->getFloat32Data();
		float *data1 = res1->getFloat32Data();

		Assert(resSize == moments.size());
		Assert(resSize == moments[0].size());

		for (int i = 0; i < resSize; i++) {
			for (int j = 0; j < resSize; j++) {
				int idx = resSize * i + j;
				data0[3 * idx + 0] = moments[i][j][0] + m_offset;
				data0[3 * idx + 1] = moments[i][j][1] + m_offset;
				
				data1[3 * idx + 0] = moments[i][j][2] + m_offset;
				data1[3 * idx + 1] = moments[i][j][3] + m_offset;
				data1[3 * idx + 2] = moments[i][j][4] + m_offset;
			}
		}

		Log(EInfo, "Output scale %d...", scale);

		char filename[256];
		sprintf(filename, "moments0_%dx.exr", scale);
		ref<FileStream> stream = new FileStream(fs::path(filename), FileStream::ETruncWrite);
		res0->write(Bitmap::EOpenEXR, stream);

		sprintf(filename, "moments1_%dx.exr", scale);
		stream = new FileStream(fs::path(filename), FileStream::ETruncWrite);
		res1->write(Bitmap::EOpenEXR, stream);

		Log(EInfo, "Output finish");
	}

	void outputBitmap(const Vector2DArray &values, char *filename) const {
		ref<Bitmap> bitmap = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, Vector2i(m_size));
		float *data = bitmap->getFloat32Data();
		for (int r = 0; r < m_size; r++) {
			for (int c = 0; c < m_size; c++) {
				*data++ = values[r][c][0];
				*data++ = values[r][c][1];
				*data++ = values[r][c][2];
			}
		}
		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		bitmap->write(Bitmap::EOpenEXR, stream);
	}

	ref<Scene> m_scene;
	Shape *m_hmap;
	ref<Sampler> m_sampler;
	int m_sqrtSpp, m_spp;
	int m_size;
	int m_maxScale;
	Float m_offset;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(LEADRplane, "Generating LEADR moment maps")
MTS_NAMESPACE_END
