#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <boost/filesystem/path.hpp>

MTS_NAMESPACE_BEGIN

struct LobeFunctor {
	LobeFunctor(char *filename, int resolution) {
		m_bitmap = new Bitmap(fs::path(filename));
		m_resolution = resolution;
	}

	Float operator()(int r, int c) const {
		int colIdx = r * m_resolution + c;
		return m_bitmap->getPixel(Point2i(colIdx, m_row))[m_channel];
	}

	void setRow(int row) {
		m_row = row;
	}

	void setChannel(int channel) {
		m_channel = channel;
	}

	ref<Bitmap> m_bitmap;
	int m_resolution;
	int m_row;
	int m_channel;
};

struct VectorFunctor {
	VectorFunctor(std::vector<Float> &vec, int resolution) {
		m_data = vec.data();
		m_resolution = resolution;
	}

	Float operator()(int r, int c) const {
		return m_data[r * m_resolution + c];
	}

	Float *m_data;
	int m_resolution;
};

class SHProjection4D : public Utility {
public:
	int run(int argc, char **argv) {
		m_bands = std::atoi(argv[1]);
		m_resolution = std::atoi(argv[2]);
		LobeFunctor fLobe(argv[3], m_resolution);
		char *filename = argv[4];
		char *reconFilename = argv[5];
		
		fLobe.setChannel(0);
		int numCoeffs = m_bands * m_bands;
		m_sh2d.resize(numCoeffs);
		for (int i = 0; i < numCoeffs; i++) {
			m_sh2d[i].clear();
		}

		// first pass
		for (int i = 0; i < m_resolution * m_resolution; i++) {
			SHVector shLobe(m_bands);
			fLobe.setRow(i);
			shLobe.projectLobe(fLobe, m_resolution);

			int k = 0;
			for (int l = 0; l < m_bands; l++) {
				for (int m = -l; m <= l; m++) {
					m_sh2d[k].push_back(shLobe(l, m));
					k++;
				}
			}
		}

		// second pass
		m_shMat.resize(numCoeffs);
		for (int i = 0; i < numCoeffs; i++) {
			m_shMat[i].resize(numCoeffs);
		}
		ref<Bitmap> shMat = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, Vector2i(numCoeffs));
		float *shData = shMat->getFloat32Data();
		for (int k1 = 0; k1 < numCoeffs; k1++) {
			VectorFunctor vec(m_sh2d[k1], m_resolution);
			SHVector shLobe(m_bands);
			shLobe.projectLobe(vec, m_resolution);

			int k2 = 0;
			for (int l = 0; l < m_bands; l++) {
				for (int m = -l; m <= l; m++) {
					float value = shLobe(l, m);
					m_shMat[k1][k2] = value;

					int idx = (k1 * numCoeffs + k2) * 3;
					shData[idx] = value;
					shData[idx + 1] = value;
					shData[idx + 2] = value;
					k2++;
				}
			}
		}

		Log(EInfo, "SH projection finished...");

		ref<FileStream> stream = new FileStream(fs::path(filename), FileStream::ETruncWrite);
		shMat->write(Bitmap::EOpenEXR, stream);

		recon4dFunction(reconFilename);
		Log(EInfo, "Reconstruction from SH finished...");

		return 0;
	}

	void recon4dFunction(char *filename) {
		int numCoeffs = m_bands * m_bands;
		
		// reverse second pass
		std::vector<std::vector<Float> > sh2d(numCoeffs, std::vector<Float>(m_resolution * m_resolution));
		for (int k1 = 0; k1 < numCoeffs; k1++) {
			SHVector shVec(m_bands);
			int k2 = 0;
			for (int l = 0; l < m_bands; l++) {
				for (int m = -l; m <= l; m++) {
					shVec(l, m) = m_shMat[k1][k2++];
				}
			}

			for (int r = 0; r < m_resolution; r++) {
				Float y = ((m_resolution - r - 1) + 0.5) / (double)m_resolution * 2.0 - 1.0;
				for (int c = 0; c < m_resolution; c++) {
					Float x = (c + 0.5) / (double)m_resolution * 2.0 - 1.0;
					int idx = r * m_resolution + c;
					
					if (x * x + y * y >= 1.0) {
						sh2d[k1][idx] = 0.0;
						continue;
					}
					Float z = sqrt(1.0 - x * x - y * y);
					Float value = shVec.eval(Vector(x, y, z));
					sh2d[k1][idx] = value;
				}
			}
		}

		// reserve first pass
		std::vector<std::vector<Float> > mat(m_resolution * m_resolution, std::vector<Float>(m_resolution * m_resolution));
		for (int i = 0; i < m_resolution * m_resolution; i++) {
			SHVector shVec(m_bands);
			int k = 0;
			for (int l = 0; l < m_bands; l++) {
				for (int m = -l; m <= l; m++) {
					shVec(l, m) = sh2d[k++][i];
				}
			}

			for (int r = 0; r < m_resolution; r++) {
				Float y = ((m_resolution - r - 1) + 0.5) / (double)m_resolution * 2.0 - 1.0;
				for (int c = 0; c < m_resolution; c++) {
					Float x = (c + 0.5) / (double)m_resolution * 2.0 - 1.0;
					int idx = r * m_resolution + c;

					if (x * x + y * y >= 1.0) {
						mat[i][idx] = 0.0;
						continue;
					}
					Float z = sqrt(1.0 - x * x - y * y);
					Float value = shVec.eval(Vector(x, y, z));
					mat[i][idx] = value;
				}
			}
		}

		// output
		ref<Bitmap> img = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, Vector2i(m_resolution * m_resolution));
		float *data = img->getFloat32Data();
		for (int i = 0; i < m_resolution * m_resolution; i++) {
			for (int j = 0; j < m_resolution * m_resolution; j++) {
				for (int c = 0; c < 3; c++) {
					int idx = ((i * m_resolution * m_resolution) + j) * 3 + c;
					data[idx] = std::max(0.0, mat[i][j]);
				}
			}
		}
		ref<FileStream> stream = new FileStream(fs::path(filename), FileStream::ETruncWrite);
		img->write(Bitmap::EOpenEXR, stream);
	}

	int m_bands;
	int m_resolution;

	// first pass: 4D -> 2D SH coefficients
	// (numCoeffs * (res^2))
	std::vector<std::vector<Float> > m_sh2d;
	// second pass: 2D -> 0D
	// (numCoeffs * numCoeffs)
	std::vector<std::vector<Float> > m_shMat;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(SHProjection4D, "Project 4D spherical function to SH coefficient matrix")
MTS_NAMESPACE_END