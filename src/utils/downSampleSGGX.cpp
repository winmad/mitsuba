/*
	Add by Lifan Wu
	Mar 20, 2016
*/

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/util.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include "../phase/microflake_fiber.h"

MTS_NAMESPACE_BEGIN

class DownSampleSGGX : public Utility {
public:
	typedef std::vector<std::vector<std::vector<Vector> > > GridData;

	typedef TPoint3<double> Pt;
	typedef TVector3<double> Vec;

	int run(int argc, char **argv) {
		if(argc != 9) {
			cout << "Down-sample SGGX volume data by a scale" << endl;
			cout << "Syntax: mtsutil downSampleVolume <orientation_volume> <density_volume> <scale> <stddev> <num_SGGX_lobes> <s1_prefix> <s2_prefix> <cdf_prefix>" << endl;
			return -1;
		}

		char *end_ptr = NULL;
		int scaleValue = strtol(argv[3], &end_ptr, 10);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse integer value");
		scale = Vector3i(scaleValue);

		Float stddev = strtod(argv[4], &end_ptr);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse floating point value");
		d = GaussianFiberDistribution(stddev);

		numSGGXlobes = strtol(argv[5], &end_ptr, 10);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse integer value");

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

		bbox = originVol->getAABB();

		/*
		Properties props2("gridvolume");
		props2.setString("filename", argv[2]);
		props2.setBoolean("sendData", false);

		VolumeDataSource *densityVol = static_cast<VolumeDataSource *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(VolumeDataSource), props2));
		densityVol->configure();
		*/

		init(originVol);
		cluster(originVol);

		bool lazy = true;
		downSample(originVol, lazy);

		Log(EInfo, "finish down-sampling, save volume data to file");
		
		for (int i = 0; i < numSGGXlobes; i++) {
			std::string s1Filename(argv[6]);
			std::string s2Filename(argv[7]);
			std::string cdfFilename(argv[8]);

			s1Filename += formatString("_%i.vol", i);
			s2Filename += formatString("_%i.vol", i);
			cdfFilename += formatString("_%i.vol", i);

			ref<FileStream> outFile = new FileStream(s1Filename.c_str(), FileStream::ETruncReadWrite);
			writeVolume(s1[i], bbox, 3, outFile);

			outFile = new FileStream(s2Filename.c_str(), FileStream::ETruncReadWrite);
			writeVolume(s2[i], bbox, 3, outFile);

			outFile = new FileStream(cdfFilename.c_str(), FileStream::ETruncReadWrite);
			writeVolume(cdf[i], bbox, 1, outFile);
		}
		
		return 0;
	}

	void init(VolumeDataSource *ori) {
		s1.resize(numSGGXlobes);
		s2.resize(numSGGXlobes);
		cdf.resize(numSGGXlobes);

		originRes = ori->getResolution();

		res.x = originRes.x / scale.x;
		res.y = originRes.y / scale.y;
		res.z = originRes.z / scale.z;

		originRes.x -= originRes.x % scale.x;
		originRes.y -= originRes.y % scale.y;
		originRes.z -= originRes.z % scale.z;

		for (int i = 0; i < numSGGXlobes; i++) {
			initS(s1[i], res);
			initS(s2[i], res);
			initS(cdf[i], res);
		}
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

	Vector6 getSValue(const Vector &w) {
		if (w.isZero()) {
			return Vector6(0.f);
		}

		Float gauss_sigma3 = d.sigmaT(1.f) * 2.f;
		Float gauss_sigma1 = d.sigmaT(0.f) * 2.f;

		Frame dFrame(w);
		Vector w1 = dFrame.s;
		Vector w2 = dFrame.t;

		Float sigma1 = gauss_sigma1;
		Float sigma2 = gauss_sigma1;
		Float sigma3 = gauss_sigma3;

		Matrix3x3 basis(w1, w2, w);
		Matrix3x3 D(Vector(sigma1 * sigma1, 0, 0), Vector(0, sigma2 * sigma2, 0), Vector(0, 0, sigma3 * sigma3));
		Matrix3x3 basisT;
		basis.transpose(basisT);
		Matrix3x3 S = basis * D * basisT;

		Vector6 res;
		res[0] = S.m[0][0]; res[1] = S.m[1][1]; res[2] = S.m[2][2];
		res[3] = S.m[0][1]; res[4] = S.m[0][2]; res[5] = S.m[1][2];
		return res;
	}

	void cluster(VolumeDataSource *ori) {
		N = scale.x * scale.y * scale.z;
		data.resize(N);
		idx.resize(N);
		segs.resize(N);

		clusterRes.resize(originRes.x);
		for (int i = 0; i < originRes.x; i++) {
			clusterRes[i].resize(originRes.y);
			for (int j = 0; j < originRes.y; j++) {
				clusterRes[i][j].resize(originRes.z);
				for (int k = 0; k < originRes.z; k++) {
					clusterRes[i][j][k] = -1;
				}
			}
		}

		Log(EInfo, "Start Clustering...");

		for (int i = 0; i < originRes.x; i += scale.x) {
			for (int j = 0; j < originRes.y; j += scale.y) {
				for (int k = 0; k < originRes.z; k += scale.z) {
					numColors = 0;
					hasht.clear();
					diffIndices.clear();
					numClusters = numSGGXlobes;

					int count = 0;

					for (int dx = 0; dx < scale.x; dx++) {
						for (int dy = 0; dy < scale.y; dy++) {
							for (int dz = 0; dz < scale.z; dz++) {
								Pt p;
								for (int c = 0; c < 3; c++) {
									p[c] = ori->lookupFloat(i + dx, j + dy, k + dz, c);
								}

								if (p.x * p.x + p.y * p.y + p.z * p.z < 1e-6)
									continue;

								int maxComponent = 0;
								double maxValue = abs(p[0]);
								for (int c = 1; c < 3; c++) {
									if (maxValue < abs(p[c])) {
										maxValue = abs(p[c]);
										maxComponent = c;
									}
								}

								if (p[maxComponent] < 0)
									p = -p;

								data[count] = p;
								idx[count] = Vector3i(i + dx, j + dy, k + dz);
								segs[count] = -1;

								int hashValue = hashFunc(p);
								if (numColors < 200 && hasht.find(hashValue) == hasht.end()) {
									hasht[hashValue] = numColors;
									diffIndices.push_back(count);
									numColors++;
								}

								count++;
							}
						}
					}

					N = count;
					numClusters = std::min(numClusters, numColors);

					//Log(EInfo, "Data prepared, total data = %d, maximum cluster = %d...", N, numClusters);

					initCluster();
					kMeans();

					//Log(EInfo, "Finish one clustering...");

					for (int c = 0; c < N; c++) {
						Vector3i nowIdx = idx[c];
						clusterRes[nowIdx.x][nowIdx.y][nowIdx.z] = segs[c];
					}
				}
			}
		}

		Log(EInfo, "Finish Clustering...\nStart Downsampling SGGX...");
		/*
		GridData tmp;
		initS(tmp, originRes);
		for (int i = 0; i < originRes.x; i++) {
			for (int j = 0; j < originRes.y; j++) {
				for (int k = 0; k < originRes.z; k++) {
					if (clusterRes[i][j][k] < 0)
						tmp[i][j][k] = Vector(0.f);
					else if (clusterRes[i][j][k] == 0)
						tmp[i][j][k] = Vector(1.f, 0.f, 0.f);
					else
						tmp[i][j][k] = Vector(0.f, 1.f, 0.f);
				}
			}
		}
		ref<FileStream> outFile = new FileStream("orientation_cluster.vol", FileStream::ETruncReadWrite);
		writeVolume(tmp, bbox, 3, outFile);
		*/
	}

	void downSample(VolumeDataSource *ori, bool lazy = false) {
		std::vector<Vector6> sumS;
		std::vector<Float> count;
		std::vector<Float> pdf;
		sumS.resize(numSGGXlobes);
		count.resize(numSGGXlobes);
		pdf.resize(numSGGXlobes);
		
		Float norm = 1.f / (scale.x * scale.y * scale.z);

		for (int i = 0; i < originRes.x; i += scale.x) {
			for (int j = 0; j < originRes.y; j += scale.y) {
				for (int k = 0; k < originRes.z; k += scale.z) {
					for (int c = 0; c < numSGGXlobes; c++) {
						sumS[c] = Vector6(0.f);
						count[c] = 0.f;
						pdf[c] = 0.f;
					}
					Float nonzero = 0.f;
					Float zeros = 0.f;

					for (int dx = 0; dx < scale.x; dx++) {
						for (int dy = 0; dy < scale.y; dy++) {
							for (int dz = 0; dz < scale.z; dz++) {
								int clusterIdx = clusterRes[i + dx][j + dy][k + dz];
								if (clusterIdx < 0)
									continue;

								Vector w;
								for (int c = 0; c < 3; c++) {
									w[c] = ori->lookupFloat(i + dx, j + dy, k + dz, c);
								}
								//Float density = densityVol->lookupFloat(i + dx, j + dy, k + dz, 0);

								Vector6 s;
								s = getSValue(w);

								sumS[clusterIdx] += s;
								count[clusterIdx] += 1.f;
								nonzero += 1.f;
							}
						}
					}
					
					if (nonzero > 0.f)
						norm = 1.f / nonzero;
					else
						norm = 0.f;
					/*
					zeros = (scale.x * scale.y * scale.z) - nonzero;

					Float nonzeroLobes = 0.f;
					for (int c = 0; c < numSGGXlobes; c++) {
						if (count[c] > 0.f)
							nonzeroLobes += 1.f;
					}
					*/
					for (int c = 0; c < numSGGXlobes; c++) {
						if (count[c] > 0.f) {
							//sumS[c] /= (count[c] + zeros / nonzeroLobes);
							sumS[c] /= count[c];
						}
						pdf[c] = count[c] * norm;

						if (!lazy) {
							s1[c][i / scale.x][j / scale.y][k / scale.z] = Vector(sumS[c][0], sumS[c][1], sumS[c][2]);
							s2[c][i / scale.x][j / scale.y][k / scale.z] = Vector(sumS[c][3], sumS[c][4], sumS[c][5]);
						}
						else {
							Matrix3x3 Q;
							Float eig[3];

							Matrix3x3 S(sumS[c][0], sumS[c][3], sumS[c][4], 
								sumS[c][3], sumS[c][1], sumS[c][5], 
								sumS[c][4], sumS[c][5], sumS[c][2]);
							S.symEig(Q, eig);
							// eig[0] < eig[1] == eig[2]
							Vector w3(Q.m[0][0], Q.m[1][0], Q.m[2][0]);
							if (!w3.isZero()) {
								s1[c][i / scale.x][j / scale.y][k / scale.z] = normalize(w3);
								s2[c][i / scale.x][j / scale.y][k / scale.z] = Vector(eig[1], eig[2], eig[0]);
							}
							else {
								Log(EInfo, "No! zero orientation vector!");
							}
						}

						if (c == 0)
							cdf[c][i / scale.x][j / scale.y][k / scale.z] = Vector(pdf[c]);
						else
							cdf[c][i / scale.x][j / scale.y][k / scale.z] = cdf[c - 1][i / scale.x][j / scale.y][k / scale.z] + 
								Vector(pdf[c]);
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

	inline int hashFunc(const Pt &u) {
		return (int)(u.x * 1e2) * 10000 + (int)(u.y * 1e2) * 100 + (int)(u.z * 1e2);
	}

	inline double dist2(const Pt &u, const Pt &v) {
		return (u - v).lengthSquared();
	}

	void assignCluster(int now) {
#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			double minDist2 = 1e10;
			int k = -1;
			for (int j = 0; j < numClusters; j++) {
				double d = dist2(data[i], centers[now][j]);
				if (d < minDist2) {
					minDist2 = d;
					k = j;
				}
			}
			segs[i] = k;
		}
	}

	void updateCenters(int now) {
		for (int i = 0; i < numClusters; i++) {
			cnt[i] = 0;
			centers[1 - now][i] = Pt(0.0);
		}

		for (int i = 0; i < N; i++) {
			int k = segs[i];
			centers[1 - now][k] += data[i];
			cnt[k] += 1;
		}

		for (int i = 0; i < numClusters; i++) {
			if (cnt[i] > 0)
				centers[1 - now][i] /= (double)cnt[i];
			else
				centers[1 - now][i] = Pt(1e8);
		}

		/*
		for (int i = 0; i < numClusters; i++) {
			Log(EInfo, "cluster %d: %d", i, cnt[i]);
		}
		Log(EInfo, "total: %d", N);
		*/
	}

	void initCluster() {
		centers[0].resize(numClusters);
		centers[1].resize(numClusters);
		cnt.resize(numClusters);
		for (int i = 0; i < numClusters; i++) {
			int k;
			while (1) {
				int u = rand() % numColors;
				k = diffIndices[u];

				bool flag = true;
				for (int j = 0; j < i; j++) {
					if (dist2(centers[0][j], data[k]) < 1e-6) {
						flag = false;
						break;
					}
				}
				if (flag)
					break;
			}
			centers[0][i] = data[k];
		}
		/*
		for (int i = 0; i < numClusters; i++) {
			Log(EInfo, "init center %d, (%.6f, %.6f, %.6f)", i, centers[0][i].x, centers[0][i].y, centers[0][i].z);
		}
		*/
	}

	void kMeans() {
		int now = 0;
		int iter = 0;
		while (1) {
			assignCluster(now);
			updateCenters(now);

			// converged?
			bool flag = true;
			for (int i = 0; i < numClusters; i++) {
				if (dist2(centers[1 - now][i], centers[now][i]) > 1e-8) {
					flag = false;
					break;
				}
			}

			iter++;
			now = 1 - now;
			//Log(EInfo, "====== Iter %d ======", iter);
			for (int i = 0; i < numClusters; i++) {
				//Log(EInfo, "old center %d, (%.8f, %.8f, %.8f)", i, centers[1 - now][i].x, centers[1 - now][i].y, centers[1 - now][i].z);
				//Log(EInfo, "center %d, (%.8f, %.8f, %.8f)", i, centers[now][i].x, centers[now][i].y, centers[now][i].z);
			}

			if (iter > 100 || flag)
				break;
		}
	}

	Vector3i scale;
	GaussianFiberDistribution d;
	int numSGGXlobes;

	Vector3i originRes;
	Vector3i res;
	AABB bbox;

	std::vector<GridData> s1;
	std::vector<GridData> s2;
	std::vector<GridData> cdf;

	std::vector<std::vector<std::vector<int> > > clusterRes;

	int N, numClusters;
	std::vector<int> segs;
	std::vector<Vector3i> idx;
	std::vector<int> cnt;
	std::vector<Pt> centers[2];
	std::vector<Pt> data;

	int numColors;
	std::map<int, int> hasht;
	std::vector<int> diffIndices;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(DownSampleSGGX, "Down-sample SGGX phase function")
MTS_NAMESPACE_END