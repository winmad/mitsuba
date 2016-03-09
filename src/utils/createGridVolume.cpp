/*
	Add by Lifan Wu
	Dec 08, 2015
*/

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/util.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

class CreateGridVolume : public Utility {
public:
	typedef std::vector<std::vector<std::vector<Vector> > > GridData;

	inline Float gaussian(Float sqrSigma, int dx, int dy, int dz) const {
		Float sqrDist = dx * dx + dy * dy + dz * dz;
		return std::exp(-sqrDist / (2 * sqrSigma));
	}

	int run(int argc, char **argv) {
		if (argc != 3) {
			cout << "Create a new grid volume" << endl;
			cout << "Syntax: mtsutil createGridVolume <existing_grid_volume_1> <existing_grid_volume_2> <new_grid_volume>" << endl;
			return -1;
		}

		Properties props("gridvolume");
		props.setString("filename", argv[1]);
		props.setBoolean("sendData", false);

		VolumeDataSource *originVol = static_cast<VolumeDataSource *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(VolumeDataSource), props));
		originVol->configure();

		AABB bbox = originVol->getAABB();
		Vector3i res = originVol->getResolution();
		int channels = 3;

		Log(EInfo, "%s", originVol->getClass()->getName().c_str());
		Log(EInfo, "res = (%d, %d, %d)", originVol->getResolution().x, originVol->getResolution().y, originVol->getResolution().z);
		Log(EInfo, "bbox = (%.6f, %.6f, %.6f), (%.6f, %.6f, %.6f)", bbox.min.x, bbox.min.y, bbox.min.z,
			bbox.max.x, bbox.max.y, bbox.max.z);
		/*
		Properties props2("gridvolume");
		props2.setString("filename", argv[2]);
		props2.setBoolean("sendData", false);

		VolumeDataSource *vol2 = static_cast<VolumeDataSource *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(VolumeDataSource), props2));
		vol2->configure();
		*/
		 
		// creating albedo

		Vector left(0.98, 0.1, 0.1);
		Vector middle(0.95, 0.1, 0.1);
		Vector right(0.8, 0.1, 0.1);

		GridData s;
		initS(s, res);

#pragma omp parallel for
		for (int i = 0; i < res.x; i++) {
			for (int j = 0; j < res.y; j++) {
				for (int k = 0; k < res.z; k++) {
					Float x = (Float)(i) / (Float)(res.x - 1.f);

					if (x < 0.5f)
						s[i][j][k] = left;
					else
						s[i][j][k] = right;

					/*
					if (x < 0.5f) {
						Float u = x * 2;
						s[i][j][k] = left * (1.f - u * u) + middle * u * u;
					}
					else {
						Float u = (x - 0.5f) * 2 - 1.f;
						s[i][j][k] = middle * u * u + right * (1.f - u * u);
					}
					*/
// 					x = 1.f - x;
// 					if (x < 0.2f)
// 						s[i][j][k] = Vector(0.99f, 0.3f, 0.3f);
// 					else if (x < 0.4f)
// 						s[i][j][k] = Vector(0.3f, 0.99f, 0.3f);
// 					else if (x < 0.6f)
// 						s[i][j][k] = Vector(0.6f, 0.1f, 0.9f);
// 					else if (x < 0.8f)
// 						s[i][j][k] = Vector(0.9f, 0.45f, 0.1f);
// 					else
// 						s[i][j][k] = Vector(0.2f, 0.2f, 0.99f);
				}
			}
		}
		

		/*
		// density filtering

		GridData density;
		initS(density, res);
#pragma omp parallel for
		for (int i = 0; i < res.x; i++) {
			for (int j = 0; j < res.y; j++) {
				for (int k = 0; k < res.z; k++) {
					Float v = originVol->lookupFloat(i, j, k, 0);
					density[i][j][k] = Vector(v);
				}
			}
		}

		GridData s;
		initS(s, res);

		Float sigma = 1.f;
		Float sqrSigma = sigma * sigma;
#pragma omp parallel for
		for (int i = 0; i < res.x; i++) {
			for (int j = 0; j < res.y; j++) {
				for (int k = 0; k < res.z; k++) {
					Float ratio = (res.x - 1.f - i) / (Float)(res.x - 1.f);
					int filterSize;
					
					//filterSize = math::roundToInt(ratio * 0.f + (1.f - ratio) * 4.f);
					if (ratio > 0.6)
						filterSize = 0;
					else
						filterSize = 4;

					Float sumDensity = 0.f;
					Float sumWeight = 0.f;
					for (int dx = -filterSize; dx <= filterSize; dx++) {
						for (int dy = -filterSize; dy <= filterSize; dy++) {
							for (int dz = -filterSize; dz <= filterSize; dz++) {
								if (i + dx < 0 || i + dx >= res.x ||
									j + dy < 0 || j + dy >= res.y ||
									k + dz < 0 || k + dz >= res.z)
									continue;
								Float weight = 1.f;//gaussian(sqrSigma, dx, dy, dz);
								sumDensity += density[i + dx][j + dy][k + dz][0] * weight;
								sumWeight += weight;
							}
						}
					}

					if (sumWeight > 0.f)
						sumDensity /= sumWeight;
					else
						sumDensity = 0.f;

					s[i][j][k] = Vector(sumDensity);
				}
			}
		}
		channels = 1;
		*/

		/*
		channels = 1;
#pragma omp parallel for
		for (int i = 0; i < res.x; i++) {
			for (int j = 0; j < res.y; j++) {
				for (int k = 0; k < res.z; k++) {
					Vector v(originVol->lookupFloat(i, j, k, 0),
						originVol->lookupFloat(i, j, k, 1),
						originVol->lookupFloat(i, j, k, 2));

					if (v.x == 0 && v.y == 0 && v.z == 0) {
						s[i][j][k] = Vector(3.f);
					}
					else if (v.z > v.x && v.x > v.y) {
						s[i][j][k] = Vector(0.f);
					}
					else if (v.x > v.y && v.y > v.z) {
						s[i][j][k] = Vector(1.f);
					}
					else {
						s[i][j][k] = Vector(2.f);
					}
				}
			}
		}
		*/
		/*
		channels = 1;
#pragma omp parallel for
		for (int i = 0; i < res.x; i++) {
			for (int j = 0; j < res.y; j++) {
				for (int k = 0; k < res.z; k++) {
					Vector v(originVol->lookupFloat(i, j, k, 0),
						originVol->lookupFloat(i, j, k, 1),
						originVol->lookupFloat(i, j, k, 2));

					if (v.x == 0 && v.y == 0 && v.z == 0) {
						s[i][j][k] = Vector(3.f);
					}
					else if (v.z > v.x && v.x > v.y) {
						s[i][j][k] = Vector(0.f);
					}
					else if (v.x > v.y && v.y > v.z) {
						s[i][j][k] = Vector(1.f);
					}
					else {
						s[i][j][k] = Vector(2.f);
					}
				}
			}
		}
		*/


		/*
#pragma omp parallel for
		for (int i = 0; i < res.x; i++) {
			for (int j = 0; j < res.y; j++) {
				for (int k = 0; k < res.z; k++) {
					Vector v(originVol->lookupFloat(i, j, k, 0),
						originVol->lookupFloat(i, j, k, 1),
						originVol->lookupFloat(i, j, k, 2));

					if (fabs(v.x) > fabs(v.y) & fabs(v.x) > fabs(v.z)) {
						if (v.x > 0.f)
							s[i][j][k] = Vector(1.f, 0.f, 0.f);
						else
							s[i][j][k] = Vector(-1.f, 0.f, 0.f);
					}
					else if (fabs(v.y) > fabs(v.x) && fabs(v.y) > fabs(v.z)) {
						if (v.y > 0)
							s[i][j][k] = Vector(0.f, 1.f, 0.f);
						else
							s[i][j][k] = Vector(0.f, -1.f, 0.f);
					}
					else {
						s[i][j][k] = Vector(0.f);
					}
				}
			}
		}
		*/

		/*
		GridData s;
		initS(s, Vector3i(2, 2, 2));
		s[1][1][1] = s[0][1][0] = s[0][0][1] = s[1][0][0] = Vector(0.95, 0.1, 0.1);
		s[0][0][0] = s[0][1][1] = s[1][1][0] = s[1][0][1] = Vector(0.95, 0.64, 0.37);

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				for (int k = 0; k < 2; k++) {
					s[i][j][k] += Vector(rand() % 10000 * 1e-5f, rand() % 10000 * 1e-5f, rand() % 10000 * 1e-5f);
					int index = (k * 2 + j) * 2 + i;
					Log(EInfo, "index %d: (%.6f, %.6f, %.6f)", index, s[i][j][k].x, s[i][j][k].y, s[i][j][k].z);
				}
		*/

// 		GridData albedo;
// 		initS(albedo, Vector3i(2, 2, 1));
		
		// fov: 1, scale: 2
// 		albedo[0][0][0] = Vector(0.935257f, 0.09f, 0.09f);
// 		albedo[0][1][0] = Vector(0.935776f, 0.08855f, 0.088547f);
// 		albedo[1][0][0] = Vector(0.937f, 0.09212f, 0.09212f);
// 		albedo[1][1][0] = Vector(0.9356f, 0.087784f, 0.08778f);

		// fov: 1, scale: 4
// 		albedo[0][0][0] = Vector(0.915451f, 0.078598f, 0.078598f);
// 		albedo[0][1][0] = Vector(0.9172516f, 0.07698f, 0.07698f);
// 		albedo[1][0][0] = Vector(0.91661f, 0.07933f, 0.07933f);
// 		albedo[1][1][0] = Vector(0.917435f, 0.077835f, 0.077835f);

		// fov: 4, scale: 2
// 		albedo[0][0][0] = Vector(0.93547f, 0.090195f, 0.090195f);
// 		albedo[0][1][0] = Vector(0.93605f, 0.088486f, 0.088486f);
// 		albedo[1][0][0] = Vector(0.93724f, 0.092f, 0.092f);
// 		albedo[1][1][0] = Vector(0.935f, 0.087387f, 0.087387f);

		// fov: 4, scale: 4
// 		albedo[0][0][0] = Vector(0.91591f, 0.079f, 0.079f);
// 		albedo[0][1][0] = Vector(0.917f, 0.07634f, 0.07634f);
// 		albedo[1][0][0] = Vector(0.91673f, 0.0797f, 0.0797f);
// 		albedo[1][1][0] = Vector(0.9167f, 0.07655f, 0.07655f);
// 
// 		Vector3i newRes = res;
// 		GridData s;
// 		initS(s, newRes);
// 
// 		AABB newBBox = bbox;
// 
// 		for (int i = 0; i < 2; i++) {
// 			for (int j = 0; j < 2; j++) {
// 				int imin = i * res[0] / 2; // 0 or res.x/2
// 				int imax = (i + 1) * res[0] / 2 - i; // res.x/2 or res.x-1
// 				int jmin = j * res[1] / 2; // 0 or res.y/2
// 				int jmax = (j + 1) * res[1] / 2 - j; // res.y/2 or res.y-1
// 
// 				for (int x = imin; x <= imax; x++) {
// 					for (int y = jmin; y <= jmax; y++) {
// 						for (int z = 0; z < newRes.z; z++) {
// 							s[x][y][z] = albedo[i][j][0];
// 						}
// 					}
// 				}
// 			}
// 		}

// 		GridData s;
// 		initS(s, res);
// 
// 		for (int i = 0; i < res[0]; i++) {
// 			for (int j = 0; j < res[1]; j++) {
// 				for (int k = 0; k < res[2]; k++) {
// 					Vector v(originVol->lookupFloat(i, j, k, 0),
// 						originVol->lookupFloat(i, j, k, 1),
// 						originVol->lookupFloat(i, j, k, 2));
// 					s[i][j][k] = v * 0.8f;
// 				}
// 			}
// 		}

		ref<FileStream> outFile = new FileStream(argv[argc - 1], FileStream::ETruncReadWrite);
		//writeVolume(s, newBBox, channels, outFile);
		writeVolume(s, bbox, channels, outFile);

		return 0;
	}

	void initS(GridData &s, const Vector3i &res) {
		s.resize(res.x);
		for (int i = 0; i < res.x; i++) {
			s[i].resize(res.y);
			for (int j = 0; j < res.y; j++) {
				s[i][j].resize(res.z);
			}
		}
	}

	void writeVolume(GridData &s, AABB &bbox, int channels, Stream *stream) {
		stream->writeChar('V');
		stream->writeChar('O');
		stream->writeChar('L');
		stream->writeChar(3);

		int type = 1;
		stream->writeInt(type);

		stream->writeInt(s.size());
		stream->writeInt(s[0].size());
		stream->writeInt(s[0][0].size());

		stream->writeInt(channels);

		stream->writeSingle(bbox.min.x);
		stream->writeSingle(bbox.min.y);
		stream->writeSingle(bbox.min.z);
		stream->writeSingle(bbox.max.x);
		stream->writeSingle(bbox.max.y);
		stream->writeSingle(bbox.max.z);

		float *data = new float[s.size() * s[0].size() * s[0][0].size() * channels];
		for (int z = 0; z < s[0][0].size(); z++) {
			for (int y = 0; y < s[0].size(); y++) {
				for (int x = 0; x < s.size(); x++) {
					for (int c = 0; c < channels; c++) {
						data[((z * s[0].size() + y) * s.size() + x) * channels + c] = s[x][y][z][c];
					}
				}
			}
		}
		stream->writeSingleArray(data, s.size() * s[0].size() * s[0][0].size() * channels);
		delete[] data;
	}

	Float m_threshold;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(CreateGridVolume, "Create grid volume")
MTS_NAMESPACE_END