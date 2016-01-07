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

class AlbedoRemapping : public Utility {
public:
	typedef Point Pt;
	typedef Vector Vec;
	typedef std::vector<std::vector<std::vector<Vector> > > GridData;

	int run(int argc, char **argv) {
		if (argc != 4) {
			cout << "Remapping albedo values according to segmentation" << endl;
			cout << "Syntax: mtsutil albedoRemapping <albedo_volume> <segmentation_volume> <new_albedo_volume>" << endl;
			return -1;
		}

		Properties props("gridvolume");
		props.setString("filename", argv[1]);
		props.setBoolean("sendData", false);

		VolumeDataSource *albedoVol = static_cast<VolumeDataSource *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(VolumeDataSource), props));
		albedoVol->configure();

		Properties props1("gridvolume");
		props1.setString("filename", argv[2]);
		props1.setBoolean("sendData", false);

		VolumeDataSource *segVol = static_cast<VolumeDataSource *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(VolumeDataSource), props1));
		segVol->configure();

		Log(EInfo, "%s", albedoVol->getClass()->getName().c_str());
		Log(EInfo, "res = (%d, %d, %d)", albedoVol->getResolution().x, albedoVol->getResolution().y, albedoVol->getResolution().z);
		Log(EInfo, "channels = %d", albedoVol->getChannels());
		Log(EInfo, "min = (%.6f, %.6f, %.6f)", albedoVol->getAABB().min.x, albedoVol->getAABB().min.y, albedoVol->getAABB().min.z);
		Log(EInfo, "max = (%.6f, %.6f, %.6f)", albedoVol->getAABB().max.x, albedoVol->getAABB().max.y, albedoVol->getAABB().max.z);

		AABB bbox = albedoVol->getAABB();
		res = albedoVol->getResolution();

		initS(s, res);

		for (int i = 0; i < res.x; i++) {
			for (int j = 0; j < res.y; j++) {
				for (int k = 0; k < res.z; k++) {
					Vector color(albedoVol->lookupFloat(i, j, k, 0),
						albedoVol->lookupFloat(i, j, k, 1),
						albedoVol->lookupFloat(i, j, k, 2));
					Float clusterIndex = segVol->lookupFloat(i, j, k, 0);

					if (segVol->lookupFloat(j, i, k, 0) == 1) {
						s[i][j][k] = Vector(0.7, 0.7, 0.975143);
					}
					else if (clusterIndex == 0) {
						s[i][j][k] = Vector(0.98453, 0.99315, 0.974124);
					} 
					else {
						s[i][j][k] = Vector(0.7, 0.952381, 0.7);
					}
				}
			}
		}

		ref<FileStream> outFile = new FileStream(argv[3], FileStream::ETruncReadWrite);
		writeVolume(bbox, 3, outFile);

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

	void writeVolume(AABB &bbox, int channels, Stream *stream) {
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

	int numClusters;
	Vector3i res;
	
	GridData s;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(AlbedoRemapping, "Remapping albedo values according to segmentation")
MTS_NAMESPACE_END