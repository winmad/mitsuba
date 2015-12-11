/*
	Add by Lifan Wu
	Dec 07, 2015
*/

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/util.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

class SplitGridVolume : public Utility {
public:
	typedef std::vector<std::vector<std::vector<Vector> > > GridData;

	int run(int argc, char **argv) {
		if (argc != 3) {
			cout << "Split grid volume data into 2*2*1 sub-volumes" << endl;
			cout << "Syntax: mtsutil splitGridVolume <grid_volume> <target_volume>" << endl;
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
		int channels = originVol->getChannels();

		Log(EInfo, "%s", originVol->getClass()->getName().c_str());
		Log(EInfo, "res = (%d, %d, %d)", originVol->getResolution().x, originVol->getResolution().y, originVol->getResolution().z);
		Log(EInfo, "channels = %d", originVol->getChannels());
		Log(EInfo, "bbox = (%.3f, %.3f, %.3f), (%.3f, %.3f, %.3f)", bbox.min.x, bbox.min.y, bbox.min.z,
			bbox.max.x, bbox.max.y, bbox.max.z);

		Vector step = bbox.getExtents();
		step[0] /= (res[0] - 1.f); step[1] /= (res[1] - 1.f); step[2] /= (res[2] - 1.f);
	
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				int imin = i * res[0] / 2; // 0 or res.x/2
				int imax = (i + 1) * res[0] / 2 - i; // res.x/2 or res.x-1
				int jmin = j * res[1] / 2; // 0 or res.y/2
				int jmax = (j + 1) * res[1] / 2 - j; // res.y/2 or res.y-1

				AABB newBBox;
				newBBox.min = Point(bbox.min.x + step.x * imin, bbox.min.y + step.y * jmin, bbox.min.z);
				newBBox.max = Point(bbox.min.x + step.x * imax, bbox.min.y + step.y * jmax, bbox.max.z);

				GridData s;
				Vector3i newRes;
				if (res[0] % 2 == 1)
					newRes.x = res[0] / 2 + 1;
				else
					newRes.x = res[0] / 2 - i + 1;
				if (res[1] % 2 == 1)
					newRes.y = res[1] / 2 + 1;
				else
					newRes.y = res[1] / 2 - j + 1;
				newRes.z = res.z;

				initS(s, newRes);

				for (int x = imin; x <= imax; x++)
				{
					for (int y = jmin; y <= jmax; y++)
					{
						for (int z = 0; z < newRes.z; z++)
						{
							if (channels == 1) {
								float v = originVol->lookupFloat(x, y, z, 0);
								s[x - imin][y - jmin][z] = Vector(v);
							}
							else {
								Vector v(originVol->lookupFloat(x, y, z, 0),
									originVol->lookupFloat(x, y, z, 1),
									originVol->lookupFloat(x, y, z, 2));
								s[x - imin][y - jmin][z] = v;
							}

						}
					}
				}

				Log(EInfo, "===== sub-volume %d =====", i * 2 + j);
				Log(EInfo, "res = (%d, %d, %d)", newRes.x, newRes.y, newRes.z);
				Log(EInfo, "bbox = (%.3f, %.3f, %.3f), (%.3f, %.3f, %.3f)", newBBox.min.x, newBBox.min.y, newBBox.min.z,
					newBBox.max.x, newBBox.max.y, newBBox.max.z);

				std::string outputFilename(argv[2]);
				outputFilename.push_back('_');
				outputFilename.push_back('0' + (i * 2 + j));
				outputFilename += ".vol";

				ref<FileStream> outFile = new FileStream(outputFilename.c_str(), FileStream::ETruncReadWrite);
				writeVolume(s, newBBox, channels, outFile);
			}
		}
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

MTS_EXPORT_UTILITY(SplitGridVolume, "Split grid volume")
MTS_NAMESPACE_END