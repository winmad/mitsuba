/*
	Add by Lifan Wu
	Nov 8, 2015
*/

#include <mitsuba/core/matrix.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/phase.h>
#include "../phase/microflake_fiber.h"

MTS_NAMESPACE_BEGIN

class microflake2SGGX : public Utility {
public:
	typedef std::vector<std::vector<std::vector<Vector> > > SData;

	int run(int argc, char **argv) {
		if (argc != 6) {
			cout << "Convert microflake phase function to SGGX phase function" << endl;
			cout << "Syntax: mtsutil microflake2SGGX 0 <x> <y> <z> <stddev>" << endl;
			cout << "Syntax: mtsutil microflake2SGGX 1 <grid_orientation_volume> <stddev> <S1_volume> <S2_volume>" << endl;
			cout << "Syntax: mtsutil microflake2SGGX 2 <hgrid_orientation_volume> <stddev> <S1_volume> <S2_volume>" << endl;
			return -1;
		}

		if (strcmp(argv[1], "0") == 0) {
			char *end_ptr = NULL;
			w3.x = strtod(argv[2], &end_ptr);
			if (*end_ptr != '\0')
				SLog(EError, "Could not parse floating point value");
			w3.y = strtod(argv[3], &end_ptr);
			if (*end_ptr != '\0')
				SLog(EError, "Could not parse floating point value");
			w3.z = strtod(argv[4], &end_ptr);
			if (*end_ptr != '\0')
				SLog(EError, "Could not parse floating point value");
			stddev = strtod(argv[5], &end_ptr);
			if (*end_ptr != '\0')
				SLog(EError, "Could not parse floating point value");

			GaussianFiberDistribution d(stddev);

			w3 = normalize(w3);
			Frame dFrame(w3);
			Vector w1 = dFrame.s;
			Vector w2 = dFrame.t;
			Log(EInfo, "Frame = ");
			Log(EInfo, "%.6f %.6f %.6f", w3.x, w3.y, w3.z);
			Log(EInfo, "%.6f %.6f %.6f", w1.x, w1.y, w1.z);
			Log(EInfo, "%.6f %.6f %.6f", w2.x, w2.y, w2.z);

			Float normFactor;
			Float sigma3 = sigma(w3, dFrame, d);
			Float sigma1 = sigma(w1, dFrame, d);
			Float sigma2 = sigma1;
// 			normFactor = std::max(sigma1, sigma3);
// 			sigma1 /= normFactor;
// 			sigma2 = sigma1;
// 			sigma3 /= normFactor;

			Float gauss_sigma3 = d.sigmaT(1.f) * 2;
			Float gauss_sigma1 = d.sigmaT(0.f) * 2;
// 			normFactor = std::max(gauss_sigma1, gauss_sigma3);
// 			gauss_sigma1 /= normFactor;
// 			gauss_sigma3 /= normFactor;

			sigma1 = sigma2 = gauss_sigma1;
			sigma3 = gauss_sigma3;

			Matrix3x3 basis(w1, w2, w3);
			Matrix3x3 D(Vector(sigma1 * sigma1, 0, 0), Vector(0, sigma2 * sigma2, 0), Vector(0, 0, sigma3 * sigma3));
			Matrix3x3 basisT;
			basis.transpose(basisT);
			Matrix3x3 S;
			S = basis * D * basisT;

			Log(EInfo, "Gaussian sigma, sigma_3 = %.6f, sigma_1 = sigma_2 = %.6f", gauss_sigma3, gauss_sigma1);
			Log(EInfo, "sigma_3 = %.6f, sigma_1 = sigma_2 = %.6f", sigma3, sigma1);
			Log(EInfo, "S = ");
			Log(EInfo, "%.6f %.6f %.6f", S.m[0][0], S.m[0][1], S.m[0][2]);
			Log(EInfo, "%.6f %.6f %.6f", S.m[1][0], S.m[1][1], S.m[1][2]);
			Log(EInfo, "%.6f %.6f %.6f", S.m[2][0], S.m[2][1], S.m[2][2]);

			Log(EInfo, "===== test property of sigma =====");

			int numTheta = 180;
			int numPhi = 360;

			Float dTheta = (M_PI / (Float)numTheta);
			Float dPhi = (M_PI * 2.f / (Float)numPhi);

			Float result = 0.f;

// 			Float sigmaT_SGGX = sigma(w3, S.m[0][0], S.m[1][1], S.m[2][2], S.m[0][1], S.m[0][2], S.m[1][2]);
// 			Float sigmaT_microflake = 2.f * d.sigmaT(1.f);
// 			Log(EInfo, "SGGX = %.6f, microflake = %.6f", sigmaT_SGGX, sigmaT_microflake);
// 			result += fabsf(sigmaT_microflake - sigmaT_SGGX);
			
			Float maxSigmaT = 0.f;

			for (int i = 0; i <= numTheta; i++) {
				Float theta = dTheta * i;
				Float z = cosf(theta);
				Float sinTheta = sinf(theta);
				for (int j = 0; j < (i == 0 ? 1 : 360); j++) {
					Float phi = dPhi * j;
					Float x = sinTheta * cosf(phi);
					Float y = sinTheta * sinf(phi);
					Vector wi = Vector(x, y, z);

					Float sigmaT_SGGX = sigma(wi, S.m[0][0], S.m[1][1], S.m[2][2], S.m[0][1], S.m[0][2], S.m[1][2]);

					Float cosTheta = dot(wi, w3);
					Float sigmaT_microflake = 2.f * d.sigmaT(cosTheta);

					maxSigmaT = std::max(maxSigmaT, sigmaT_SGGX);

					result += fabsf(sigmaT_microflake - sigmaT_SGGX);
				}
			}
			result /= numTheta * numPhi;
			Log(EInfo, "max sigmaT = %.6f", maxSigmaT);
			Log(EInfo, "average error in sigma is %.6f", result);
		}
		else if (strcmp(argv[1], "1") == 0) {
			Properties props("gridvolume");
			props.setString("filename", argv[2]);
			props.setBoolean("sendData", false);
			VolumeDataSource *ori = static_cast<VolumeDataSource *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(VolumeDataSource), props));
			ori->configure();

			char *end_ptr = NULL;
			stddev = strtod(argv[3], &end_ptr);
			if (*end_ptr != '\0')
				SLog(EError, "Could not parse floating point value");

			GaussianFiberDistribution d(stddev);

			Log(EInfo, "Start S matrix calculation");
			SData s1, s2;
			gridOri2S(ori, d, s1, s2);

			Log(EInfo, "Calculation finished. Start writing volume file");
			AABB bbox = ori->getAABB();
			ref<FileStream> s1File = new FileStream(argv[4], FileStream::ETruncReadWrite);
			ref<FileStream> s2File = new FileStream(argv[5], FileStream::ETruncReadWrite);
			gridWriteS(s1, bbox, s1File.get());
			gridWriteS(s2, bbox, s2File.get());
		}
		else if (strcmp(argv[2], "1") == 0) {

		}

		return 0;
	}

	void initS(SData &s, const Vector3i &res) {
		s.resize(res.x);
		for (int i = 0; i < res.x; i++) {
			s[i].resize(res.y);
			for (int j = 0; j < res.y; j++) {
				s[i][j].resize(res.z);
			}
		}
	}

	void getSValue(Vector &w, const GaussianFiberDistribution &d, Vector &s1Val, Vector &s2Val) {
		if (w.isZero()) {
			s1Val = Vector(0.f);
			s2Val = Vector(0.f);
			return;
		}

		Float gauss_sigma3 = d.sigmaT(1.f) * 2;
		Float gauss_sigma1 = d.sigmaT(0.f) * 2;

		Frame dFrame(w);
		Vector w1 = dFrame.s;
		Vector w2 = dFrame.t;

// 		Float sigma3 = sigma(w, dFrame, d);
// 		Float sigma1 = sigma(w1, dFrame, d);
// 		Float sigma2 = sigma1;
		
		Float sigma1 = gauss_sigma1;
		Float sigma2 = gauss_sigma1;
		Float sigma3 = gauss_sigma3;

		Matrix3x3 basis(w1, w2, w);
		Matrix3x3 D(Vector(sigma1 * sigma1, 0, 0), Vector(0, sigma2 * sigma2, 0), Vector(0, 0, sigma3 * sigma3));
		Matrix3x3 basisT;
		basis.transpose(basisT);
		Matrix3x3 S = basis * D * basisT;

		s1Val[0] = S.m[0][0]; s1Val[1] = S.m[1][1]; s1Val[2] = S.m[2][2];
		s2Val[0] = S.m[0][1]; s2Val[1] = S.m[0][2]; s2Val[2] = S.m[1][2];
	}

	void gridOri2S(VolumeDataSource *ori, const GaussianFiberDistribution &d, SData &s1, SData &s2) {
		Vector3i res = ori->getResolution();

		initS(s1, res);
		initS(s2, res);

#pragma omp parallel for
		for (int i = 0; i < res.x; i++) {
			for (int j = 0; j < res.y; j++) {
				for (int k = 0; k < res.z; k++) {
					Vector w(ori->lookupFloat(i, j, k, 0), ori->lookupFloat(i, j, k, 1), ori->lookupFloat(i, j, k, 2));
					Vector s1Val(0.f), s2Val(0.f);
					getSValue(w, d, s1Val, s2Val);
					s1[i][j][k] = s1Val;
					s2[i][j][k] = s2Val;
				}
			}
		}
	}

	void gridWriteS(SData &s, AABB &bbox, Stream *stream) {
		stream->writeChar('V');
		stream->writeChar('O');
		stream->writeChar('L');
		stream->writeChar(3);

		int type = 1;
		stream->writeInt(type);

		stream->writeInt(s.size());
		stream->writeInt(s[0].size());
		stream->writeInt(s[0][0].size());

		int channels = 3;
		stream->writeInt(channels);

		stream->writeSingle(bbox.min.x);
		stream->writeSingle(bbox.min.y);
		stream->writeSingle(bbox.min.z);
		stream->writeSingle(bbox.max.x);
		stream->writeSingle(bbox.max.y);
		stream->writeSingle(bbox.max.z);

		float *data = new float[s.size() * s[0].size() * s[0][0].size() * 3];
		for (int z = 0; z < s[0][0].size(); z++) {
			for (int y = 0; y < s[0].size(); y++) {
				for (int x = 0; x < s.size(); x++) {
					for (int c = 0; c < 3; c++) {
						data[((z * s[0].size() + y) * s.size() + x) * 3 + c] = s[x][y][z][c];
					}
				}
			}
		}
		stream->writeSingleArray(data, s.size() * s[0].size() * s[0][0].size() * 3);
		delete[] data;
	}

	// numerical integration has problems!
	Float sigma(const Vector &w, const Frame &dFrame, const GaussianFiberDistribution &d) {
		Frame frame(w);
		Float result = 0.f;

		int numTheta = 90;
		int numPhi = 360;

		Float dTheta = (M_PI * 0.5f / (Float)numTheta);
		Float dPhi = (M_PI * 2.f / (Float)numPhi);

		for (int i = 0; i <= numTheta; i++) {
			Float theta = dTheta * i;
			Float z = cosf(theta);
			Float sinTheta = sinf(theta);
			for (int j = 0; j < (i == 0 ? 1 : 360); j++) {
				Float phi = dPhi * j;
				Float x = sinTheta * cosf(phi);
				Float y = sinTheta * sinf(phi);
				
				Vector localWm = Vector(x, y, z);
				Vector wm = frame.toWorld(localWm);
				result += z * d.pdf(dFrame.toLocal(wm)) * sinTheta;
			}
		}
		result *= dTheta * dPhi;
		return result;
	}

	Float sigma(const Vector &wi, Float Sxx, Float Syy, Float Szz,
		Float Sxy, Float Sxz, Float Syz) const {
		Float sigmaSqr = wi.x * wi.x * Sxx + wi.y * wi.y * Syy + wi.z * wi.z * Szz +
			2.f * (wi.x * wi.y * Sxy + wi.x * wi.z * Sxz + wi.y * wi.z * Syz);
		return (sigmaSqr > 0.f) ? sqrtf(sigmaSqr) : 0.f;
	}

	Vector w3;
	Float stddev;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(microflake2SGGX, "Convert microflake phase function to SGGX phase function")
MTS_NAMESPACE_END