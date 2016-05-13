/*
	Add by Lifan Wu
	May 03, 2016
*/

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/util.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include "../phase/microflake_fiber.h"

MTS_NAMESPACE_BEGIN

class AlignLobes : public Utility {
public:
	typedef std::vector<std::vector<std::vector<Vector> > > GridData;

	int run(int argc, char **argv) {
		if (argc != 8) {
			cout << "Align SGGX lobes" << endl;
			cout << "Syntax: mtsutil alignLobes <num_SGGX_lobes> <s1_prefix> <s2_prefix> <cdf_prefix> <newS1_prefix> <newS2_prefix> <newCdf_prefix>" << endl;
			return -1;
		}

		char *end_ptr = NULL;
		numSGGXlobes = strtol(argv[1], &end_ptr, 10);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse integer value");

		s1.resize(numSGGXlobes);
		s2.resize(numSGGXlobes);
		cdf.resize(numSGGXlobes);
		pdf.resize(numSGGXlobes);

		bool lazy = true;

		VolumeDataSource *originVol;

		for (int i = 0; i < numSGGXlobes; i++) {
			std::string s1Filename(argv[2]);
			std::string s2Filename(argv[3]);
			std::string cdfFilename(argv[4]);

			s1Filename += formatString("_%i.vol", i);
			s2Filename += formatString("_%i.vol", i);
			cdfFilename += formatString("_%i.vol", i);

			Properties props("gridvolume");
			props.setString("filename", s1Filename);
			props.setBoolean("sendData", false);
			originVol = static_cast<VolumeDataSource *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(VolumeDataSource), props));
			originVol->configure();
			copyFromVol(s1[i], originVol);

			Properties props2("gridvolume");
			props2.setString("filename", s2Filename);
			props2.setBoolean("sendData", false);
			originVol = static_cast<VolumeDataSource *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(VolumeDataSource), props2));
			originVol->configure();
			copyFromVol(s2[i], originVol);

			Properties props3("gridvolume");
			props3.setString("filename", cdfFilename);
			props3.setBoolean("sendData", false);
			originVol = static_cast<VolumeDataSource *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(VolumeDataSource), props3));
			originVol->configure();
			copyFromVol(cdf[i], originVol);
		}

		res = originVol->getResolution();

		for (int i = 0; i < numSGGXlobes; i++) {
			initS(pdf[i], res);
			if (i == 0)
				pdf[i] = cdf[i];
			else {
				for (int x = 0; x < res.x; x++)
					for (int y = 0; y < res.y; y++)
						for (int z = 0; z < res.z; z++)
							pdf[i][x][y][z] = cdf[i][x][y][z] - cdf[i - 1][x][y][z];
			}
		}

		Log(EInfo, "%s", originVol->getClass()->getName().c_str());
		Log(EInfo, "res = (%d, %d, %d)", originVol->getResolution().x, originVol->getResolution().y, originVol->getResolution().z);
		Log(EInfo, "channels = %d", originVol->getChannels());
		Log(EInfo, "min = (%.6f, %.6f, %.6f)", originVol->getAABB().min.x, originVol->getAABB().min.y, originVol->getAABB().min.z);
		Log(EInfo, "max = (%.6f, %.6f, %.6f)", originVol->getAABB().max.x, originVol->getAABB().max.y, originVol->getAABB().max.z);

		bbox = originVol->getAABB();

		alignLobes();

		Log(EInfo, "finish aligning, save volume data to file");

		for (int i = 0; i < numSGGXlobes; i++) {
			std::string s1Filename(argv[5]);
			std::string s2Filename(argv[6]);
			std::string cdfFilename(argv[7]);

			s1Filename += formatString("_%i.vol", i);
			s2Filename += formatString("_%i.vol", i);
			cdfFilename += formatString("_%i.vol", i);

			ref<FileStream> outFile = new FileStream(s1Filename.c_str(), FileStream::ETruncReadWrite);
			writeVolume(newS1[i], bbox, 3, outFile);

			outFile = new FileStream(s2Filename.c_str(), FileStream::ETruncReadWrite);
			writeVolume(newS2[i], bbox, 3, outFile);

			outFile = new FileStream(cdfFilename.c_str(), FileStream::ETruncReadWrite);
			writeVolume(newCdf[i], bbox, 1, outFile);
		}

		return 0;
	}

	void copyFromVol(GridData &s, VolumeDataSource *ori) {
		int channels = ori->getChannels();
		res = ori->getResolution();
		initS(s, res);

		for (int i = 0; i < res.x; i++)
			for (int j = 0; j < res.y; j++)
				for (int k = 0; k < res.z; k++) {
					if (channels == 1) {
						float v = ori->lookupFloat(i, j, k, 0);
						s[i][j][k] = Vector(v);
					}
					else {
						Vector v(ori->lookupFloat(i, j, k, 0),
							ori->lookupFloat(i, j, k, 1),
							ori->lookupFloat(i, j, k, 2));
						s[i][j][k] = v;
					}
				}
	}

	void initS(GridData &s, const Vector3i &res) {
		s.resize(res.x);
		for (int i = 0; i < res.x; i++) {
			s[i].resize(res.y);
			for (int j = 0; j < res.y; j++) {
				s[i][j].resize(res.z);
				for (int k = 0; k < res.z; k++) {
					s[i][j][k] = Vector(0.f);
				}
			}
		}
	}

	void alignLobes() {
		stdVecs.clear();
		if (numSGGXlobes == 2) {
			stdVecs.push_back(Vector(1.f, 0.f, 0.f));
			stdVecs.push_back(Vector(0.f, 1.f, 0.f));
		}
		else if (numSGGXlobes == 3) {
			stdVecs.push_back(Vector(1.f, 0.f, 0.f));
			stdVecs.push_back(Vector(0.f, 1.f, 0.f));
			stdVecs.push_back(Vector(0.f, 0.f, 1.f));
		}

		newS1.resize(numSGGXlobes);
		newS2.resize(numSGGXlobes);
		newCdf.resize(numSGGXlobes);
		for (int l = 0; l < numSGGXlobes; l++) {
			initS(newS1[l], res);
			initS(newS2[l], res);
			initS(newCdf[l], res);
		}

		std::vector<bool> flag(numSGGXlobes, false);

		for (int i = 0; i < res.x; i++) {
			for (int j = 0; j < res.y; j++) {
				for (int k = 0; k < res.z; k++) {
					for (int l = 0; l < numSGGXlobes; l++)
						flag[l] = false;
					for (int l = 0; l < numSGGXlobes; l++) {
						int choose = -1;
						float maxv = 0.f;
						if (pdf[l][i][j][k].x > 0) {
							for (int c = 0; c < stdVecs.size(); c++) {
								if (flag[c])
									continue;
								float v = abs(dot(s1[l][i][j][k], stdVecs[c]));
								if (v > maxv) {
									maxv = v;
									choose = c;
								}
							}

							if (choose != -1) {
								newCdf[choose][i][j][k] = Vector(pdf[l][i][j][k]);
								newS1[choose][i][j][k] = s1[l][i][j][k];
								newS2[choose][i][j][k] = s2[l][i][j][k];
								flag[choose] = true;
							}
						}
					}

					for (int l = 1; l < numSGGXlobes; l++) {
						newCdf[l][i][j][k] += newCdf[l - 1][i][j][k];
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

	int numSGGXlobes;

	Vector3i originRes;
	Vector3i res;
	AABB bbox;

	std::vector<GridData> s1;
	std::vector<GridData> s2;
	std::vector<GridData> cdf;
	std::vector<GridData> pdf;

	std::vector<GridData> newS1;
	std::vector<GridData> newS2;
	std::vector<GridData> newCdf;

	std::vector<Vector> stdVecs;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(AlignLobes, "Align SGGX lobes")
MTS_NAMESPACE_END