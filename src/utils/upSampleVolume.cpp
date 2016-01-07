/*
	Add by Lifan Wu
	Jan 6, 2016
*/

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/util.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

class UpSampleVolume : public Utility {
public:
	typedef std::vector<std::vector<std::vector<Vector> > > GridData;

	int run(int argc, char **argv) {
		if (argc != 4) {
			cout << "Up-sample grid volume data by a scale" << endl;
			cout << "Syntax: mtsutil upSampleVolume <grid_volume> <scale> <target_volume>" << endl;
			return -1;
		}

		int scaleValue = 1;
		char *end_ptr = NULL;
		scaleValue = strtol(argv[2], &end_ptr, 10);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse integer value");
		scale = Vector3i(scaleValue);

		Properties props("gridvolume");
		props.setString("filename", argv[1]);
		props.setBoolean("sendData", false);

		VolumeDataSource *originVol = static_cast<VolumeDataSource *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(VolumeDataSource), props));
		originVol->configure();

		Log(EInfo, "%s", originVol->getClass()->getName().c_str());
		Log(EInfo, "res = (%d, %d, %d)", originVol->getResolution().x, originVol->getResolution().y, originVol->getResolution().z);
		Log(EInfo, "channels = %d", originVol->getChannels());
		Log(EInfo, "min = (%.6f, %.6f, %.6f)", originVol->getAABB().min.x, originVol->getAABB().min.y, originVol->getAABB().min.z);
		Log(EInfo, "max = (%.6f, %.6f, %.6f)", originVol->getAABB().max.x, originVol->getAABB().max.y, originVol->getAABB().max.z);

		AABB bbox = originVol->getAABB();

		GridData s;
		upSample(originVol, s);

		Log(EInfo, "finish up-sampling, save volume data to file");
		ref<FileStream> outFile = new FileStream(argv[3], FileStream::ETruncReadWrite);
		writeVolume(s, bbox, originVol->getChannels(), outFile);

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

	void upSample(VolumeDataSource *ori, GridData &s) {
		int channels = ori->getChannels();
		Vector3i originRes = ori->getResolution();

		Vector3i res;
		res.x = originRes.x * scale.x;
		res.y = originRes.y * scale.y;
		res.z = originRes.z * scale.z;

		initS(s, res);

		for (int i = 0; i < res.x; i++) {
			for (int j = 0; j < res.y; j++) {
				for (int k = 0; k < res.z; k++) {
					if (channels == 1) {
						s[i][j][k] = Vector(ori->lookupFloat(i / scale.x, j / scale.y, k / scale.z, 0));
					}
					else {
						s[i][j][k] = Vector(ori->lookupFloat(i / scale.x, j / scale.y, k / scale.z, 0),
							ori->lookupFloat(i / scale.x, j / scale.y, k / scale.z, 1),
							ori->lookupFloat(i / scale.x, j / scale.y, k / scale.z, 2));
					}
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

MTS_EXPORT_UTILITY(UpSampleVolume, "Up-sample volume")
MTS_NAMESPACE_END