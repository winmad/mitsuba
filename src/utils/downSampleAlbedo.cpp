/*
	Add by Lifan Wu
	Jan 3, 2016
*/

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/util.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

class DownSampleAlbedo : public Utility {
public:
	typedef std::vector<std::vector<std::vector<Vector> > > GridData;

	int run(int argc, char **argv) {
		if (argc != 5) {
			cout << "Down-sample albedo volume data by a scale" << endl;
			cout << "Syntax: mtsutil downSampleVolume <albedo_volume> <density_volume> <scale> <new_albedo_volume>" << endl;
			return -1;
		}

		char *end_ptr = NULL;
		int scaleValue = strtol(argv[3], &end_ptr, 10);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse integer value");
		scale = Vector3i(scaleValue);

		Properties props("gridvolume");
		props.setString("filename", argv[1]);
		props.setBoolean("sendData", false);

		VolumeDataSource *albedoVol = static_cast<VolumeDataSource *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(VolumeDataSource), props));
		albedoVol->configure();

		Properties props1("gridvolume");
		props1.setString("filename", argv[2]);
		props1.setBoolean("sendData", false);

		VolumeDataSource *densityVol = static_cast<VolumeDataSource *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(VolumeDataSource), props1));
		densityVol->configure();

		Log(EInfo, "%s", albedoVol->getClass()->getName().c_str());
		Log(EInfo, "res = (%d, %d, %d)", albedoVol->getResolution().x, albedoVol->getResolution().y, albedoVol->getResolution().z);
		Log(EInfo, "channels = %d", albedoVol->getChannels());
		Log(EInfo, "min = (%.6f, %.6f, %.6f)", albedoVol->getAABB().min.x, albedoVol->getAABB().min.y, albedoVol->getAABB().min.z);
		Log(EInfo, "max = (%.6f, %.6f, %.6f)", albedoVol->getAABB().max.x, albedoVol->getAABB().max.y, albedoVol->getAABB().max.z);

		AABB bbox = albedoVol->getAABB();

		GridData s;
		downSampleAlbedo(albedoVol, densityVol, s);

		Log(EInfo, "finish down-sampling, save volume data to file");
		ref<FileStream> outFile = new FileStream(argv[4], FileStream::ETruncReadWrite);
		writeVolume(s, bbox, 3, outFile);

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

	void downSampleAlbedo(VolumeDataSource *albedoVol, VolumeDataSource *densityVol, GridData &s) {
		int channels = 3;
		Vector3i originRes = albedoVol->getResolution();

		Vector3i res;
		res.x = originRes.x / scale.x;
		res.y = originRes.y / scale.y;
		res.z = originRes.z / scale.z;

		initS(s, res);

		originRes.x -= originRes.x % scale.x;
		originRes.y -= originRes.y % scale.y;
		originRes.z -= originRes.z % scale.z;

#pragma omp parallel for
		for (int i = 0; i < originRes.x; i += scale.x) {
			for (int j = 0; j < originRes.y; j += scale.y) {
				for (int k = 0; k < originRes.z; k += scale.z) {
					Vector sumAlbedoWeighted(0.f);
					Vector sumAlbedo(0.f);
					Float sumDensity = 0.f;
					for (int dx = 0; dx < scale.x; dx++) {
						for (int dy = 0; dy < scale.y; dy++) {
							for (int dz = 0; dz < scale.z; dz++) {
								Vector albedo(albedoVol->lookupFloat(i + dx, j + dy, k + dz, 0),
									albedoVol->lookupFloat(i + dx, j + dy, k + dz, 1),
									albedoVol->lookupFloat(i + dx, j + dy, k + dz, 2));
								Float density = densityVol->lookupFloat(i + dx, j + dy, k + dz, 0);
								sumAlbedoWeighted += albedo * density;
								sumAlbedo += albedo;
								sumDensity += density;
							}
						}
					}
					if (sumDensity > 1e-6f)
						s[i / scale.x][j / scale.y][k / scale.z] = sumAlbedoWeighted / sumDensity;
					else
						s[i / scale.x][j / scale.y][k / scale.z] = sumAlbedo / (Float)(scale.x * scale.y * scale.z);
				}
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

	Vector3i scale;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(DownSampleAlbedo, "Down-sample albedo volume")
MTS_NAMESPACE_END