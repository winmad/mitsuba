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
								p.x += step.x * (i + dx) * step.x;
								p.y += step.y * (j + dy) * step.y;
								p.z += step.z * (k + dz) * step.z;

								if (channels == 1) {
									// Only handle 1 channel
									Float v = originVol->lookupFloat(p);
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

		/*
		int xres = stream->readInt(),
			yres = stream->readInt(),
			zres = stream->readInt();
		m_res = Vector3i(xres, yres, zres);
		m_channels = stream->readInt();
		std::string format;

		switch (type) {
		case EFloat32:
			if (m_channels != 1 && m_channels != 3)
				Log(EError, "Encountered an unsupported float32 volume data "
				"file (%i channels, only 1 and 3 are supported)",
				m_channels);
			format = "float32";
			break;
		case EFloat16:
			format = "float16";
			Log(EError, "Error: float16 volumes are not yet supported!");
		case EUInt8:
			format = "uint8";
			if (m_channels != 1 && m_channels != 3)
				Log(EError, "Encountered an unsupported uint8 volume data "
				"file (%i channels, only 1 and 3 are supported)", m_channels);
			break;
		case EQuantizedDirections:
			format = "qdir";
			if (m_channels != 3)
				Log(EError, "Encountered an unsupported quantized direction "
				"volume data file (%i channels, only 3 are supported)",
				m_channels);
			break;
		default:
			Log(EError, "Encountered a volume data file of unknown type (type=%i, channels=%i)!", type, m_channels);
		}

		m_volumeType = (EVolumeType)type;

		if (!m_dataAABB.isValid()) {
			Float xmin = stream->readSingle(),
				ymin = stream->readSingle(),
				zmin = stream->readSingle();
			Float xmax = stream->readSingle(),
				ymax = stream->readSingle(),
				zmax = stream->readSingle();
			m_dataAABB = AABB(Point(xmin, ymin, zmin), Point(xmax, ymax, zmax));
		}

		Log(EDebug, "Mapped \"%s\" into memory: %ix%ix%i (%i channels, format = %s), %s, %s",
			resolved.filename().string().c_str(), m_res.x, m_res.y, m_res.z, m_channels, format.c_str(),
			memString(m_mmap->getSize()).c_str(), m_dataAABB.toString().c_str());
		m_data = (uint8_t *)(((float *)m_mmap->getData()) + 12);
		*/
	}

	Vector3i scale;
	int channels;
	AABB bbox;

	// new resolution of down-samples grid volume data
	Vector3i res;
	
	// original step size
	Vector step;

	uint8_t *data;
	std::vector<std::vector<std::vector<Float> > > gridData;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(DownSampleVolume, "Down-sample volume")
MTS_NAMESPACE_END