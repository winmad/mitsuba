/*
	Add by Lifan Wu
	Nov 8, 2015
*/

#include <mitsuba/core/matrix.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/phase.h>
#include "../phase/microflake_fiber.h"

MTS_NAMESPACE_BEGIN

class microflake2SGGX : public Utility {
public:
	int run(int argc, char **argv) {
		if (argc != 5) {
			cout << "Convert microflake phase function to SGGX phase function" << endl;
			cout << "Syntax: mtsutil microflake2SGGX <x> <y> <z> <stddev>" << endl;
			return -1;
		}

		char *end_ptr = NULL;
		w3.x = strtod(argv[1], &end_ptr);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse floating point value");
		w3.y = strtod(argv[2], &end_ptr);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse floating point value");
		w3.z = strtod(argv[3], &end_ptr);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse floating point value");
		stddev = strtod(argv[4], &end_ptr);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse floating point value");

		GaussianFiberDistribution d(stddev);
		dFrame = Frame(w3);
		Vector w1 = dFrame.s;
		Vector w2 = dFrame.t;

		Float sigma3 = sigma(w3, d);
		Float sigma1 = sigma(w1, d);
		Float sigma2 = sigma1;

		Matrix3x3 basis(w1, w2, w3);
		Matrix3x3 D(Vector(sigma1 * sigma1, 0, 0), Vector(0, sigma2 * sigma2, 0), Vector(0, 0, sigma3 * sigma3));
		Matrix3x3 basisT;
		basis.transpose(basisT);
		Matrix3x3 S = basis * D * basisT;

		Log(EInfo, "sigma_3 = %.6f, sigma_1 = sigma_2 = %.6f", sigma3, sigma1);
		Log(EInfo, "S = ");
		Log(EInfo, "%.6f %.6f %.6f", S.m[0][0], S.m[0][1], S.m[0][2]);
		Log(EInfo, "%.6f %.6f %.6f", S.m[1][0], S.m[1][1], S.m[1][2]);
		Log(EInfo, "%.6f %.6f %.6f", S.m[2][0], S.m[2][1], S.m[2][2]);
		return 0;
	}

	Float sigma(const Vector &w, const GaussianFiberDistribution &d) {
		Frame frame(w);
		Float result = 0.f;
		for (int i = 0; i <= 90; i++) {
			Float theta = M_PI / 180.f * i;
			Float z = cosf(theta);
			Float sinTheta = sinf(theta);
			for (int j = 0; j < (i == 0 ? 1 : 360); j++) {
				Float phi = M_PI / 180.f * j;
				Float x = sinTheta * cosf(phi);
				Float y = sinTheta * sinf(phi);
				
				Vector localWm = Vector(x, y, z);
				Vector wm = frame.toWorld(localWm);
				result += z * d.pdf(dFrame.toLocal(wm));
			}
		}
		result *= (M_PI * M_PI) / (180.f * 180.f);
		return result;
	}

	Vector w3;
	Float stddev;
	Frame dFrame;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(microflake2SGGX, "Convert microflake phase function to SGGX phase function")
MTS_NAMESPACE_END