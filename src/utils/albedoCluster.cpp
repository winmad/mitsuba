/*
	Add by Lifan Wu
	Dec 30, 2015
*/

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/util.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

class AlbedoCluster : public Utility {
public:
	typedef TPoint3<double> Pt;
	typedef TVector3<double> Vec;

	inline double dist2(const Pt &u, const Pt &v) {
		return (u - v).lengthSquared();
	}

	void assignCluster(int now) {
#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			double minDist2 = 1e10;
			float k = -1.f;
			for (int j = 0; j < numClusters; j++) {
				double d = dist2(data[i], centers[now][j]);
				if (d < minDist2) {
					minDist2 = d;
					k = (float)j;
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
			int k = (int)segs[i];
			centers[1 - now][k] += data[i];
			cnt[k] += 1;
		}

		for (int i = 0; i < numClusters; i++) {
			if (cnt[i] > 0)
				centers[1 - now][i] /= (double)cnt[i];
			else
				centers[1 - now][i] = Pt(1e8);
		}

		for (int i = 0; i < numClusters; i++) {
			Log(EInfo, "cluster %d: %d", i, cnt[i]);
		}
		Log(EInfo, "total: %d", N);
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

		for (int i = 0; i < numClusters; i++) {
			Log(EInfo, "init center %d, (%.6f, %.6f, %.6f)", i, centers[0][i].x, centers[0][i].y, centers[0][i].z);
		}
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
			Log(EInfo, "====== Iter %d ======", iter);
			for (int i = 0; i < numClusters; i++) {
				//Log(EInfo, "old center %d, (%.8f, %.8f, %.8f)", i, centers[1 - now][i].x, centers[1 - now][i].y, centers[1 - now][i].z);
				Log(EInfo, "center %d, (%.8f, %.8f, %.8f)", i, centers[now][i].x, centers[now][i].y, centers[now][i].z);
			}

			if (iter > 100 || flag)
				break;
		}
	}

	int run(int argc, char **argv) {
		if (argc != 4) {
			cout << "Use k-means for albedo voxels clustering" << endl;
			cout << "Syntax: mtsutil albedoCluster <albedo_volume> <number_of_clusters> <segmentation_volume>" << endl;
			return -1;
		}

		char *end_ptr = NULL;
		numClusters = strtol(argv[2], &end_ptr, 10);
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

		AABB bbox = originVol->getAABB();
		res = originVol->getResolution();
		N = res.x * res.y * res.z;

		numColors = 0;
		hasht.clear();
		diffIndices.clear();

		double spaceTerm = 0.1;

		data.resize(N);
		for (int i = 0; i < res.x; i++) {
			for (int j = 0; j < res.y; j++) {
				for (int k = 0; k < res.z; k++) {
					int index = (k * res.y + j) * res.x + i;
					Pt p;
					for (int c = 0; c < 3; c++) {
						p[c] = originVol->lookupFloat(i, j, k, c);
					}
					/*
					Point3 pos((i + 1.f) / res.x, (j + 1.f) / res.y, (k + 1.f) / res.z);
					pos *= spaceTerm;
					p[3] = pos[0]; p[4] = pos[1]; p[5] = pos[2];
					*/
					data[index] = p;

					int hashValue = hashFunc(p);
					if (numColors < 1000 && hasht.find(hashValue) == hasht.end()) {
						//Log(EInfo, "%d, (%.6f, %.6f, %.6f)", hashValue, p.x, p.y, p.z);
						hasht[hashValue] = numColors;
						diffIndices.push_back(index);
						numColors++;
					}
				}
			}
		}

		segs = new float[N];
		for (int i = 0; i < N; i++) segs[i] = -1.f;
		
		Log(EInfo, "Start Clustering...");
		
		initCluster();
		kMeans();

		ref<FileStream> outFile = new FileStream(argv[3], FileStream::ETruncReadWrite);
		writeVolume(bbox, 1, outFile);
		
		return 0;
	}

	inline int hashFunc(const Pt &u) {
		return (int)(u.x * 1e2) * 10000 + (int)(u.y * 1e2) * 100 + (int)(u.z * 1e2);
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

		stream->writeSingleArray(segs, N);
	}

	int numClusters;
	Vector3i res;
	int N;
	std::vector<int> cnt;
	std::vector<Pt> centers[2];
	std::vector<Pt> data;
	float *segs;

	int numColors;
	std::map<int, int> hasht;
	std::vector<int> diffIndices;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(AlbedoCluster, "K-means clustering for albedo voxels")
MTS_NAMESPACE_END