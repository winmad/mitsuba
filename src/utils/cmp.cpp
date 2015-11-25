/*
	Add by Lifan Wu
	Nov 12, 2015
*/

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/util.h>

MTS_NAMESPACE_BEGIN

class Cmp : public Utility {
public:
	int run(int argc, char **argv) {
		if(argc != 4) {
			cout << "Down-sample grid volume data by a scale" << endl;
			cout << "Syntax: mtsutil cmp <result_image> <reference_image> <bool_clamp>" << endl;
			return -1;
		}

		char *end_ptr = NULL;
		int clampFlag = strtol(argv[3], &end_ptr, 10);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse integer value");

		float *resImg;
		float *refImg;
		int widthRes, heightRes;
		int widthRef, heightRef;

		loadPfm(argv[1], resImg, widthRes, heightRes, clampFlag);
		loadPfm(argv[2], refImg, widthRef, heightRef, clampFlag);

		double relRmse = relRMSE(resImg, refImg, widthRes, heightRes);
		double rmse = RMSE(resImg, refImg, widthRes, heightRes);

		//Log(EInfo, "Relative RMSE = %.6f", relRmse);
		Log(EInfo, "%s", argv[1]);
		Log(EInfo, "RMSE = %.6f", rmse);
	}

	double relRMSE(float *resImg, float *refImg, int width, int height) {
		double sum = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				for (int c = 0; c < 3; c++) {
					double d = resImg[3 * (y * width + x) + c] - refImg[3 * (y * width + x) + c];
					double refVal = refImg[3 * (y * width + x) + c];
					if (refVal > 1e-4f)
						sum += d * d / (refVal * refVal);
					else
						sum += d * d;
				}
			}
		}
		sum /= width * height * 3;
		return sqrt(sum);
	}

	double RMSE(float *resImg, float *refImg, int width, int height) {
		double sum = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				for (int c = 0; c < 3; c++) {
					double d = resImg[3 * (y * width + x) + c] - refImg[3 * (y * width + x) + c];
					sum += d * d;
				}
			}
		}
		sum /= width * height * 3;
		return sqrt(sum);
	}

	void loadPfm(char *filename, float* &data, int &width, int &height, int clampFlag) {
		FILE *fp = fopen(filename, "rb");
		float maxVal;
		fscanf(fp, "PF\n%d %d\n%f\n", &width, &height, &maxVal);
		
		data = new float[width * height * 3];
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				fread((float*)&data[3 * (y * width + x) + 0], sizeof(float), 1, fp);
				fread((float*)&data[3 * (y * width + x) + 1], sizeof(float), 1, fp);
				fread((float*)&data[3 * (y * width + x) + 2], sizeof(float), 1, fp);
			}
		}
		fclose(fp);

		if (clampFlag) {
			maxVal = 1.f;
			for (int i = 0; i < width * height * 3; i++)
				data[i] = math::clamp(data[i], 0.f, maxVal);
		}
	}

	void savePfm(char *filename, float *data, int width, int height) {
		char header[512];
		sprintf(header, "PF\n%d %d\n-1.000000\n", width, height);

		FILE *fp = fopen(filename, "wb");
		fwrite(header, strlen(header), 1, fp);
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				fwrite(&data[3 * (y * width + x) + 0], sizeof(float), 1, fp);
				fwrite(&data[3 * (y * width + x) + 1], sizeof(float), 1, fp);
				fwrite(&data[3 * (y * width + x) + 2], sizeof(float), 1, fp);
			}
		}
		fclose(fp);
	}

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(Cmp, "Compare the result image with the reference image")
MTS_NAMESPACE_END