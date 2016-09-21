/*
	Add by Lifan Wu
	Sep 7, 2015
*/

#include <mitsuba/core/matrix.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include "../phase/microflake_fiber.h"

MTS_NAMESPACE_BEGIN

class calcDownSampleError : public Utility {
public:
	typedef double Float;
	typedef std::vector<std::vector<std::vector<Vector> > > GridData;

	struct MultiLobeSGGX {
		int numSGGXlobes;
		bool lazy;

		Vector3i size;
		std::vector<GridData> s1;
		std::vector<GridData> s2;
		std::vector<GridData> cdf;

		MultiLobeSGGX(int _numSGGXlobes, bool _lazy)
			: numSGGXlobes(_numSGGXlobes), lazy(_lazy) {
			s1.resize(numSGGXlobes);
			s2.resize(numSGGXlobes);
			cdf.resize(numSGGXlobes);
			size = Vector3i(0);
		}

		void loadVolume(const std::string& filename, GridData* gridData) {
			Properties props("gridvolume");
			props.setString("filename", filename);
			props.setBoolean("sendData", false);
			VolumeDataSource* vol = static_cast<VolumeDataSource*>(PluginManager::getInstance()->
				createObject(MTS_CLASS(VolumeDataSource), props));
			vol->configure();

			if (size.x == 0)
				size = vol->getResolution();

			gridData->resize(size.x);

#pragma omp parallel for
			for (int i = 0; i < size.x; i++) {
				(*gridData)[i].resize(size.y);
				for (int j = 0; j < size.y; j++) {
					(*gridData)[i][j].resize(size.z);
					for (int k = 0; k < size.z; k++) {
						if (vol->getChannels() == 1) {
							Float value = vol->lookupFloat(i, j, k, 0);
							(*gridData)[i][j][k] = Vector3(value);
						}
						else {
							for (int c = 0; c < 3; c++) {
								(*gridData)[i][j][k][c] = vol->lookupFloat(i, j, k, c);
							}
						}
					}
				}
			}
		}

		void loadSGGX(const std::string& dataDir, int downScale) {
			for (int i = 0; i < numSGGXlobes; i++) {
				std::string filename = dataDir + formatString("cdf_%ix_%i.vol", downScale, i);
				loadVolume(filename, &cdf[i]);
				filename = dataDir + formatString("s1_%ix_%i.vol", downScale, i);
				loadVolume(filename, &s1[i]);
				filename = dataDir + formatString("s2_%ix_%i.vol", downScale, i);
				loadVolume(filename, &s2[i]);
			}

// 			Log(EInfo, "resolution = (%d, %d, %d)", size.x, size.y, size.z);
// 			int x = 8;
// 			int y = 6; 
// 			int z = 94;
// 			Log(EInfo, "s1[0] value at (%d, %d, %d): (%.6f, %.6f, %.6f)", x, y, z,
// 				s1[0][x][y][z][0],
// 				s1[0][x][y][z][1],
// 				s1[0][x][y][z][2]);
// 			Log(EInfo, "s2[0] value at (%d, %d, %d): (%.6f, %.6f, %.6f)", x, y, z,
// 				s2[0][x][y][z][0],
// 				s2[0][x][y][z][1],
// 				s2[0][x][y][z][2]);
		}

		Vector6 getSValue(const Vector &w, const Vector &sigma) {
			if (w.isZero() || sigma.isZero()) {
				return Vector6(0.f);
			}

			Frame dFrame(w);
			Vector w1 = dFrame.s;
			Vector w2 = dFrame.t;

			Float sigma1 = sigma[0];
			Float sigma2 = sigma[1];
			Float sigma3 = sigma[2];

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

		void configure() {
			if (!lazy)
				return;

#pragma omp parallel for
			for (int l = 0; l < numSGGXlobes; l++) {
				for (int i = 0; i < size.x; i++) {
					for (int j = 0; j < size.y; j++) {
						for (int k = 0; k < size.z; k++) {
							Vector w;
							Vector sigma;
							for (int c = 0; c < 3; c++) {
								w[c] = s1[l][i][j][k][c];
								sigma[c] = s2[l][i][j][k][c];
							}
							
							if (!w.isZero())
								w = normalize(w);
							Vector6 sggxValue = getSValue(w, sigma);
							
							for (int c = 0; c < 3; c++) {
								s1[l][i][j][k][c] = sggxValue[c];
								s2[l][i][j][k][c] = sggxValue[c + 3];
							}
						}
					}
				}
			}

// 			int x = 8;
// 			int y = 6;
// 			int z = 94;
// 			Log(EInfo, "s1[0] value at (%d, %d, %d): (%.6f, %.6f, %.6f)", x, y, z,
// 				s1[0][x][y][z][0],
// 				s1[0][x][y][z][1],
// 				s1[0][x][y][z][2]);
// 			Log(EInfo, "s2[0] value at (%d, %d, %d): (%.6f, %.6f, %.6f)", x, y, z,
// 				s2[0][x][y][z][0],
// 				s2[0][x][y][z][1],
// 				s2[0][x][y][z][2]);
		}

		Float ndf(const Vector &wm, Float Sxx, Float Syy, Float Szz,
			Float Sxy, Float Sxz, Float Syz) const {
			Float detS = Sxx * Syy * Szz - Sxx * Syz * Syz - Syy * Sxz * Sxz - Szz * Sxy * Sxy + 2.f * Sxy * Sxz * Syz;
			Float den = wm.x * wm.x * (Syy * Szz - Syz * Syz) + wm.y * wm.y * (Sxx * Szz - Sxz * Sxz) +
				wm.z * wm.z * (Sxx * Syy - Sxy * Sxy) + 2.f * (wm.x * wm.y * (Sxz * Syz - Szz * Sxy) +
				wm.x * wm.z * (Sxy * Syz - Syy * Sxz) + wm.y * wm.z * (Sxy * Sxz - Sxx * Syz));
			//SLog(EInfo, "[%.6f, %.6f]", detS, den);
			return powf(fabsf(detS), 1.5) / (M_PI * den * den);
		}

		Float getNdf(const Vector& w, int l, int x, int y, int z) {
			Float res = 0.f;

			Float Sxx = s1[l][x][y][z][0];
			Float Syy = s1[l][x][y][z][1];
			Float Szz = s1[l][x][y][z][2];
			Float Sxy = s2[l][x][y][z][0];
			Float Sxz = s2[l][x][y][z][1];
			Float Syz = s2[l][x][y][z][2];
				
			return ndf(w, Sxx, Syy, Szz, Sxy, Sxz, Syz);
		}

		bool isEmpty(int x, int y, int z) {
			return (cdf[numSGGXlobes - 1][x][y][z][0] < 0.5f);
		}
	};

	Float getNdfDiff(Vector3i posRef, Vector3i posDown, int downScale) {
		const int numTheta = 60;
		const int numPhi = 90;

		Float dTheta = M_PI / ((Float)numTheta + 1);
		Float dPhi = 2.f * M_PI / (Float)numPhi;

		Float refValues[numTheta + 1][numPhi];
		Float singleRef[numTheta + 1][numPhi];

		Vector3 s1, s2;

		for (int i = 0; i <= numTheta; i++) {
			for (int j = 0; j < numPhi; j++) {
				refValues[i][j] = 0.f;
			}
		}

		Float norm_factor = 1e2f;

		Float cnt = 0.f;
		for (int x = posRef.x; x < posRef.x + downScale; x++) {
			for (int y = posRef.y; y < posRef.y + downScale; y++) {
				for (int z = posRef.z; z < posRef.z + downScale; z++) {
					if (refPhase->isEmpty(x, y, z))
						continue;
					cnt += 1.f;
					
					// get single distribution
					Float sumValue = 0.f;
					for (int i = 0; i <= numTheta; i++) {
						Float theta = dTheta * i;
						Float vz = cos(theta);
						Float sinTheta = sin(theta);
						for (int j = 0; j < numPhi; j++) {
							Float phi = dPhi * j;
							Float vx = sinTheta * cos(phi);
							Float vy = sinTheta * sin(phi);
							Vector wi(vx, vy, vz);

							singleRef[i][j] = refPhase->getNdf(wi, 0, x, y, z);
							sumValue += singleRef[i][j];
						}
					}

					if (sumValue < 1e-6f)
						Log(EInfo, "bad singleRef");

					// normalize
					for (int i = 0; i <= numTheta; i++) {
						for (int j = 0; j < numPhi; j++) {
							singleRef[i][j] /= sumValue;
							singleRef[i][j] *= norm_factor;
							refValues[i][j] += singleRef[i][j];
						}
					}
				}
			}
		}

		if (cnt < 1e-6)
			Log(EInfo, "bad refValues");

		for (int i = 0; i <= numTheta; i++) {
			for (int j = 0; j < numPhi; j++) {
				refValues[i][j] /= cnt;
			}
		}

		Float downValues[numTheta + 1][numPhi];
		Float singleDown[numTheta + 1][numPhi];

		for (int i = 0; i <= numTheta; i++) {
			for (int j = 0; j < numPhi; j++) {
				downValues[i][j] = 0.f;
			}
		}

		for (int l = 0; l < numSGGXlobes; l++) {
			int x = posDown.x;
			int y = posDown.y; 
			int z = posDown.z;
			Float pdf = downPhase->cdf[l][x][y][z][0];
			if (l > 0)
				pdf -= downPhase->cdf[l - 1][x][y][z][0];

			if (pdf < 1e-6f)
				continue;

			// get single distribution
			Float sumValue = 0.f;
			for (int i = 0; i <= numTheta; i++) {
				Float theta = dTheta * i;
				Float vz = cos(theta);
				Float sinTheta = sin(theta);
				for (int j = 0; j < numPhi; j++) {
					Float phi = dPhi * j;
					Float vx = sinTheta * cos(phi);
					Float vy = sinTheta * sin(phi);
					Vector wi(vx, vy, vz);

					singleDown[i][j] = downPhase->getNdf(wi, l, x, y, z);
					sumValue += singleDown[i][j];
				}
			}

			if (sumValue < 1e-6)
				Log(EInfo, "bad singleDown");

			// normalize
			for (int i = 0; i <= numTheta; i++) {
				for (int j = 0; j < numPhi; j++) {
					singleDown[i][j] /= sumValue;
					singleDown[i][j] *= norm_factor;
					downValues[i][j] += pdf * singleDown[i][j];
				}
			}
		}

		Float res = 0.0;

		Float sumRef = 0.0;
		Float sumDown = 0.f;

		for (int i = 0; i <= numTheta; i++) {
			Float theta = dTheta * i;
			Float sinTheta = sin(theta);
			for (int j = 0; j < numPhi; j++) {
				//SLog(EInfo, "%.6f, %.6f", refValues[i][j], downValues[i][j]);
				res += pow(downValues[i][j] - refValues[i][j], 2.0) * sinTheta;
				//res += abs(downValues[i][j] - refValues[i][j]) * sinTheta;

				sumRef += refValues[i][j];
				sumDown += downValues[i][j];
			}
		}
		res *= dTheta * dPhi;

		if (fabs(sumRef - norm_factor) > 1 || fabs(sumDown - norm_factor) > 1)
			Log(EInfo, "not normalized: %.6f, %.6f", sumRef, sumDown);

		return res;
	}

	Float getError() {
		std::vector<std::vector<std::vector<Float> > > errorValues;
		std::vector<std::vector<std::vector<bool> > > isEmpty;

		errorValues.resize(downPhase->size.x);
		isEmpty.resize(downPhase->size.x);
		for (int i = 0; i < downPhase->size.x; i++) {
			errorValues[i].resize(downPhase->size.y);
			isEmpty[i].resize(downPhase->size.y);
			for (int j = 0; j < downPhase->size.y; j++) {
				errorValues[i][j].resize(downPhase->size.z);
				isEmpty[i][j].resize(downPhase->size.z);
			}
		}

#pragma omp parallel for
		for (int i = 0; i < downPhase->size.x; i++) {
			for (int j = 0; j < downPhase->size.y; j++) {
				for (int k = 0; k < downPhase->size.z; k++) {
					if (downPhase->isEmpty(i, j, k)) {
						isEmpty[i][j][k] = true;
						continue;
					}
					
					isEmpty[i][j][k] = false;
					Vector3i posDown(i, j, k);
					Vector3i posRef(i * downScale, j * downScale, k * downScale);
					errorValues[i][j][k] = getNdfDiff(posRef, posDown, downScale);
				}
			}
		}

		Float cnt = 0.f;
		Float res = 0.f;
		for (int i = 0; i < downPhase->size.x; i++) {
			for (int j = 0; j < downPhase->size.y; j++) {
				for (int k = 0; k < downPhase->size.z; k++) {
					if (isEmpty[i][j][k])
						continue;
					res += errorValues[i][j][k];
					cnt += 1.f;
				}
			}
		}

		Log(EInfo, "%.6f / %.6f", res, cnt);
		return res / cnt;
	}

	int run(int argc, char **argv) {
		if (argc != 5) {
			cout << "Calculate the error between reference SGGX and downsampled one" << endl;
			cout << "Syntax: mtsutil calcDownSampleError <ref_sggx_dir> <down_sggx_dir> <scale> <num_sggx_lobes>";
			return -1;
		}

		char *end_ptr = NULL;
		downScale = strtol(argv[3], &end_ptr, 10);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse integer value");
		
		numSGGXlobes = strtol(argv[4], &end_ptr, 10);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse integer value");

		bool lazy = false;
		refPhase = new MultiLobeSGGX(1, lazy);
		downPhase = new MultiLobeSGGX(numSGGXlobes, lazy);

		refPhase->loadSGGX(argv[1], 1);
		refPhase->configure();
		SLog(EInfo, "Reference phase function loaded.");
		
		downPhase->loadSGGX(argv[2], downScale);
		downPhase->configure();
		SLog(EInfo, "Downsampled phase function loaded.");

		Float error = getError();
		Log(EInfo, "Downsample phase approximation error: %.6f", error);

		return 0;
	}

	int downScale;
	int numSGGXlobes;

	MultiLobeSGGX* refPhase;
	MultiLobeSGGX* downPhase;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(calcDownSampleError, "Calculate the error between reference SGGX and downsampled one")
MTS_NAMESPACE_END