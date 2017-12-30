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

class TestTextureBitmap : public Utility {
public:
	typedef std::vector<std::vector<Vector3d> > Vector2DArray;

	int run(int argc, char **argv) {
		m_size = std::atoi(argv[2]);

		Properties props;
		props = Properties("bitmap");
		props.setString("filename", argv[1]);
		props.setString("wrapMode", "repeat");
		props.setString("filterType", "nearest");
		props.setFloat("gamma", 1.0);
		props.setFloat("uscale", 1.0);
		props.setFloat("vscale", 1.0);
		ref<Texture> tex = static_cast<Texture *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Texture), props));
		tex->configure();

		Float step = 1.0 / 512.0;
		Intersection its;
		its.hasUVPartials = false;

		Vector2DArray values(m_size, std::vector<Vector3d>(m_size));

		for (int y = 0; y < m_size; y++) {
			for (int x = 0; x < m_size; x++) {
				its.uv = Point2(step * (x + 0.25), step * (y + 0.75));
				Spectrum tmp = tex->eval(its);
				for (int c = 0; c < 3; c++)
					values[y][x][c] = tmp[c];
			}
		}
		/*
		its.uv = Point2(1.65 * step, 1.35 * step);
		Spectrum res;
		res = tex->eval(its);
		Log(EInfo, "%.6f, %.6f, %.6f", res[0], res[1], res[2]);
		*/

		outputBitmap(values, argv[3]);

		return 0;
	}

	void outputBitmap(const Vector2DArray &bsdfValues, char *filename) const {
		ref<Bitmap> bitmap = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, Vector2i(m_size));
		float *data = bitmap->getFloat32Data();
		for (int r = 0; r < m_size; r++) {
			for (int c = 0; c < m_size; c++) {
				*data++ = bsdfValues[r][c][0];
				*data++ = bsdfValues[r][c][1];
				*data++ = bsdfValues[r][c][2];
			}
		}
		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		bitmap->write(Bitmap::EOpenEXR, stream);
	}

	int m_size;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(TestTextureBitmap, "Visualize BSDF eval")
MTS_NAMESPACE_END