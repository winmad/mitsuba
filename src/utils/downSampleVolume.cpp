/*
	Add by Lifan Wu
	Nov 3, 2015
*/

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/util.h>

MTS_NAMESPACE_BEGIN

class DownSampleVolume : public Utility {
public:
	int run(int argc, char **argv) {
		if(argc != 6) {
			cout << "Down-sample grid volume data by a scale" << endl;
			cout << "Syntax: mtsutil downSampleVolume <origin_volume> <scale x> <scale y> <scale z> <target_volume>" << endl;
			return -1;
		}

		char *end_ptr = NULL;
		scale.x = strtol(argv[2], &end_ptr, 10);
		scale.y = strtol(argv[3], &end_ptr, 10);
		scale.z = strtol(argv[4], &end_ptr, 10);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse integer value");

		ref<FileStream> outFile = new FileStream(argv[5], FileStream::ETruncReadWrite);

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

		channels = originVol->getChannels();
		bbox = originVol->getAABB();
		
		Vector3i originRes = originVol->getResolution();
		res.x = originRes.x / scale.x;
		res.y = originRes.y / scale.y;
		res.z = originRes.z / scale.z;

		gridData.resize(res.x);
		for (int i = 0; i < res.x; i++) {
			gridData[i].resize(res.y);
			for (int j = 0; j < res.y; j++) {
				gridData[i][j].resize(res.z);
			}
		}

		Vector extents(bbox.getExtents());
		step = Vector(extents.x / (originRes.x - 1.f), extents.y / (originRes.y - 1.f), extents.z / (originRes.z - 1.f));

		originRes.x -= originRes.x % scale.x;
		originRes.y -= originRes.y % scale.y;
		originRes.z -= originRes.z % scale.z;

		for (int i = 0; i < originRes.x; i += scale.x) {
			for (int j = 0; j < originRes.y; j += scale.y) {
				for (int k = 0; k < originRes.z; k += scale.z) {
					Float sumValue = 0.f;
					Float cnt = 0.f;
					for (int dx = 0; dx < scale.x; dx++) {
						for (int dy = 0; dy < scale.y; dy++) {
							for (int dz = 0; dz < scale.z; dz++) {
								Point p = bbox.min;
								p.x += step.x * (i + dx);
								p.y += step.y * (j + dy);
								p.z += step.z * (k + dz);

								if (channels == 1) {
									// Only handle 1 channel
									float v = originVol->lookupFloat(p);
									sumValue += v;
									cnt += 1.f;
								}
								else {
									Log(EError, "Cannot handle 3 channels...");
								}
							}
						}
					}
					gridData[i / scale.x][j / scale.y][k / scale.z] = sumValue / cnt;
				}
			}
		}

		Log(EInfo, "finish down-sampling, save volume data to file");
		saveToFile(outFile);

		return 0;
	}

	void saveToFile(Stream *stream) {
		stream->writeChar('V');
		stream->writeChar('O');
		stream->writeChar('L');
		stream->writeChar(3);
		
		int type = 1;
		stream->writeInt(type);

		stream->writeInt(res.x);
		stream->writeInt(res.y);
		stream->writeInt(res.z);

		stream->writeInt(channels);

		stream->writeSingle(bbox.min.x);
		stream->writeSingle(bbox.min.y);
		stream->writeSingle(bbox.min.z);
		stream->writeSingle(bbox.max.x);
		stream->writeSingle(bbox.max.y);
		stream->writeSingle(bbox.max.z);

		data = new float[res.x * res.y * res.z];
		for (int z = 0; z < res.z; z++) {
			for (int y = 0; y < res.y; y++) {
				for (int x = 0; x < res.x; x++) {
					data[(z * res.y + y) * res.x + x] = gridData[x][y][z];
				}
			}
		}
		stream->writeSingleArray(data, res.x * res.y * res.z);
	}

	Vector3i scale;
	int channels;
	AABB bbox;

	// new resolution of down-samples grid volume data
	Vector3i res;
	
	// original step size
	Vector step;

	float *data;
	std::vector<std::vector<std::vector<float> > > gridData;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(DownSampleVolume, "Down-sample volume")
MTS_NAMESPACE_END