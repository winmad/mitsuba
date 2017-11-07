#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/math.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/sampler.h>
#include <boost/filesystem/path.hpp>

MTS_NAMESPACE_BEGIN

class BeckmannSurface : public Utility {
public:
	int run(int argc, char **argv) {
		if (argc != 6) {
			cout << "Generate Beckmann surface" << endl;
			cout << "Syntax: mtsutil beckmannSurface <alpha_x> <alpha_y> <N> <size> <seed>" << endl;
			return -1;
		}

		char *end_ptr = NULL;
		m_alpha_x = strtod(argv[1], &end_ptr);
		m_alpha_y = strtod(argv[2], &end_ptr);
		m_N = strtol(argv[3], &end_ptr, 10);
		m_size = strtol(argv[4], &end_ptr, 10);
		m_seed = strtol(argv[5], &end_ptr, 10);

		Properties props = Properties("independent");
		props.setInteger("seed", m_seed);
		Sampler *sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), props));
		sampler->configure();

		init(sampler);
		genHeightmap();
		genNormalmap();

		return 0;
	}

	void init(Sampler *sampler) {
		m_phi.resize(m_N);
		m_fx.resize(m_N);
		m_fy.resize(m_N);

		m_scale = sqrt(2.0 / (double)m_N);
		for (int i = 0; i < m_N; i++) {
			Float U1 = sampler->next1D();
			Float U2 = sampler->next1D();
			Float U3 = sampler->next1D();

			m_phi[i] = 2.0 * M_PI * U1;
			double theta = 2.0 * M_PI * U2;
			double r = sqrt(-math::fastlog(U3));

			m_fx[i] = r * std::cos(theta) * m_alpha_x;
			m_fy[i] = r * std::sin(theta) * m_alpha_y;
		}

		m_pos.resize(m_size);
		double step = 2.0 / (m_size - 1);
		for (int i = 0; i < m_size; i++)
			m_pos[i] = -1.0 + step * i;
	}

	void genHeightmap() {
		ref<Bitmap> bitmap = new Bitmap(Bitmap::EPixelFormat::ELuminance, Bitmap::EFloat32,
			Vector2i(m_size), -1);
		float *data = bitmap->getFloat32Data();
		float min_height = 1e4f;

#pragma omp parallel for
		for (int k = 0; k < m_size * m_size; k++) {
			int c = k % m_size;
			int r = k / m_size;

			double x = c;
			double y = r;

			double res = 0.0;
			for (int i = 0; i < m_N; i++) {
				res += std::cos(m_phi[i] + x * m_fx[i] + y * m_fy[i]);
			}
			data[k] = res * m_scale;
		}

		for (int k = 0; k < m_size * m_size; k++)
			min_height = std::min(min_height, data[k]);
		for (int k = 0; k < m_size * m_size; k++)
			data[k] -= min_height;

		fs::path filename("heightmap_2d.exr");
		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		bitmap->write(Bitmap::EOpenEXR, stream);
	}

	void genNormalmap() {
		ref<Bitmap> slope = new Bitmap(Bitmap::EPixelFormat::ERGB, Bitmap::EFloat32,
			Vector2i(m_size), -1);
		ref<Bitmap> normal = new Bitmap(Bitmap::EPixelFormat::ERGB, Bitmap::EFloat32,
			Vector2i(m_size), -1);
		float *slopeData = slope->getFloat32Data();
		float *normalData = normal->getFloat32Data();

#pragma omp parallel for
		for (int k = 0; k < m_size * m_size; k++) {
			int c = k % m_size;
			int r = k / m_size;

			double x = c;
			double y = r;

			Vector3d res(0.0);
			Vector3d ndir(0.0);
			for (int i = 0; i < m_N; i++) {
				double tp = std::sin(m_phi[i] + x * m_fx[i] + y * m_fy[i]);
				res.x += -m_fx[i] * tp;
				res.y += -m_fy[i] * tp;
			}
			res *= m_scale;
			
			ndir = Vector3d(-res.x, -res.y, 1.0);
			ndir = normalize(ndir);

			for (int i = 0; i < 3; i++) {
				slopeData[k * 3 + i] = res[i];
				normalData[k * 3 + i] = ndir[i];
			}
		}

		fs::path filename("normalmap_2d.exr");
		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		normal->write(Bitmap::EOpenEXR, stream);

		filename = fs::path("slope_2d.exr");
		stream = new FileStream(filename, FileStream::ETruncWrite);
		slope->write(Bitmap::EOpenEXR, stream);
	}

	double m_alpha_x, m_alpha_y;
	double m_scale;
	int m_N, m_size;
	uint64_t m_seed;

	std::vector<double> m_phi, m_fx, m_fy;
	std::vector<double> m_pos;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(BeckmannSurface, "Beckmann surface generator")
MTS_NAMESPACE_END
