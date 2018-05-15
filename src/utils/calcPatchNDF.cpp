#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/vmf.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <boost/filesystem/path.hpp>
#include "../bsdfs/microfacet.h"
#include "multi_vmf.h"

MTS_NAMESPACE_BEGIN

class CalcPatchNDF : public Utility {
public:
	int run(int argc, char **argv) {
		m_wi = Vector(std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3]));
		m_maskingOption = std::atoi(argv[4]);
		m_scene = loadScene(argv[5]);
		m_size = std::atoi(argv[6]);
		m_sqrtSpp = std::atoi(argv[7]);
		Float xmin = std::atof(argv[8]);
		Float xmax = std::atof(argv[9]);
		Float ymin = std::atof(argv[10]);
		Float ymax = std::atof(argv[11]);
		m_normal_stLim = std::atof(argv[12]);
		fs::path filename(argv[13]);

		m_wi = normalize(m_wi);
		m_aabb = AABB2(Point2(xmin, ymin), Point2(xmax, ymax));
		m_scene->initialize();
		m_spp = m_sqrtSpp * m_sqrtSpp;

		Properties props = Properties("independent");
		props.setInteger("seed", 19931004);
		m_samplers.resize(233);
		m_samplers[0] = static_cast<Sampler *>(PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));
		m_samplers[0]->configure();
		for (int i = 1; i < 233; i++) {
			m_samplers[i] = m_samplers[0]->clone();
		}

		// ground truth NDF
		std::vector<std::vector<Vector> > normals(m_sqrtSpp, std::vector<Vector>(m_sqrtSpp));
		std::vector<std::vector<double> > weights(m_sqrtSpp, std::vector<double>(m_sqrtSpp));
		
#pragma omp parallel for
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				normals[i][j] = sampleNormal(i, j, weights[i][j]);
				normals[i][j] /= m_normal_stLim;
			}
		}
		Log(EInfo, "Finish binning.");

		double normFactor = (double)(m_size * m_size) / (4.0 * m_spp);
		Normal avgNormal(0.0);
		double dp = 4.0 / (double)(m_size * m_size);
		double totWeight = 0.0f;
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				if (weights[i][j] > 0) {
					avgNormal += normals[i][j] / normals[i][j].z;
					totWeight += weights[i][j] / normals[i][j].z;
				}
			}
		}
		avgNormal /= (double)m_spp;
		totWeight /= (double)m_spp;	
		Log(EInfo, "Validation: %.8f should equal to %.8f", totWeight, dot(avgNormal, m_wi));

		m_D = new Bitmap(Bitmap::ELuminance, Bitmap::EFloat32, Vector2i(m_size, m_size));
		m_D->clear();
		float *dData = m_D->getFloat32Data();
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				Vector &normal = normals[i][j];
				int c = math::clamp(math::floorToInt((normal.x + 1.0) * 0.5 * m_size), 0, m_size - 1);
				int r = math::clamp(math::floorToInt((normal.y + 1.0) * 0.5 * m_size), 0, m_size - 1);

				dData[(m_size - r - 1) * m_size + c] += 1.0 * normFactor;
			}
		}

		// ensure \int D(wm)<wm,wg>dwm = 1
		double totProjArea = 0.0;
		for (int r = 0; r < m_size; r++) {
			double y = (r + 0.5) / (double)m_size * 2.0 - 1.0;
			for (int c = 0; c < m_size; c++) {
				double x = (c + 0.5) / (double)m_size * 2.0 - 1.0;
				double sinTheta2 = x * x + y * y;
				if (sinTheta2 >= 1.0)
					continue;
				double jacobian = 1.0 / std::sqrt(1.0 - sinTheta2);
				totProjArea += dData[(m_size - r - 1) * m_size + c] * dp;  //std::sqrt(1.0 - sinTheta2) * jacobian
			}
		}
		Log(EInfo, "Validation: int[D(wm)cos(wm)] = %.8f should be equal to 1", totProjArea);

		m_res = new Bitmap(Bitmap::ELuminance, Bitmap::EFloat32, Vector2i(m_size, m_size));
		m_res->clear();
		float *data = m_res->getFloat32Data();

		m_hmap = m_scene->getShapes()[0];
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				Vector &normal = normals[i][j];

				//int c = (normal.x > 0.9999f ? m_size - 1 : math::floorToInt((normal.x + 1.0) * 0.5 * m_size));
				//int r = (normal.y > 0.9999f ? m_size - 1 : math::floorToInt((normal.y + 1.0) * 0.5 * m_size));
				int c = math::clamp(math::floorToInt((normal.x + 1.0) * 0.5 * m_size), 0, m_size - 1);
				int r = math::clamp(math::floorToInt((normal.y + 1.0) * 0.5 * m_size), 0, m_size - 1);

// 				if (r == 64 && c == 64) {
// 					Log(EInfo, "(%.6f, %.6f, %.6f), (%d, %d), %.6f", normal.x, normal.y, normal.z,
// 						i, j, normFactor * (weights[i][j] / totWeight));
// 				}

				data[(m_size - r - 1) * m_size + c] += normFactor * (weights[i][j] / totWeight);
			}
		}
		
		// ensure \int D_{wi}(wm)dwm = 1
		double totD = 0.0;
		for (int r = 0; r < m_size; r++) {
			double y = (r + 0.5) / (double)m_size * 2.0 - 1.0;
			for (int c = 0; c < m_size; c++) {
				double x = (c + 0.5) / (double)m_size * 2.0 - 1.0;
				double sinTheta2 = x * x + y * y;
				if (sinTheta2 >= 1.0)
					continue;
				double jacobian = 1.0 / std::sqrt(1.0 - sinTheta2);
				totD += data[(m_size - r - 1) * m_size + c] * jacobian * dp;
			}
		}
		Log(EInfo, "Validation: int[D_wi(wm)] = %.8f should be equal to 1", totD);

		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		m_res->write(Bitmap::EOpenEXR, stream);

		// output NDF derived by LEADR
		if (false) {
			char fname[256];
			sprintf(fname, "LEADR_%s", argv[13]);
			outputLeadrNDF(normals, fname);
		}

		// fitting with multi-lobe vMF
		if (argc > 14) {
			m_numLobes = std::atoi(argv[14]);
			initKMeans(normals, weights);

			char fname[256];
			sprintf(fname, "vmf_init_%s", argv[13]);
			m_vmfs.outputDistribution(m_size, fname);

			int maxIters = 20;
			for (int i = 0; i < maxIters; i++) {
				EM(normals, weights);

 				Log(EInfo, "====== Iter %d vMF distribution ======", i);
				for (int l = 0; l < m_numLobes; l++) {
					Log(EInfo, "alpha = %.6f, kappa = %.6f, mu = (%.6f, %.6f, %.6f)",
						m_vmfs.m_alpha[l], m_vmfs.m_dist[l].getKappa(),
						m_vmfs.m_mu[l].x, m_vmfs.m_mu[l].y, m_vmfs.m_mu[l].z);
				}
			}

			sprintf(fname, "vmf_final_%s", argv[13]);
			m_vmfs.outputDistribution(m_size, fname);

			FILE *fp = fopen("vmf_params.txt", "w");
			fprintf(fp, "%d\n", m_numLobes);
			for (int l = 0; l < m_numLobes; l++) {
				fprintf(fp, "%.8f %.8f %.8f %.8f %.8f\n", 
					m_vmfs.m_alpha[l], m_vmfs.m_dist[l].getKappa(),
					m_vmfs.m_mu[l].x, m_vmfs.m_mu[l].y, m_vmfs.m_mu[l].z);
			}
			fclose(fp);
		}

		// output slopes
		/*
		char fname[256];
		sprintf(fname, "slopes_%s", argv[13]);
		int len = strlen(fname);
		fname[len - 3] = 't';
		fname[len - 2] = 'x';
		fname[len - 1] = 't';
		FILE *fp = fopen(fname, "w");
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				if (weights[i][j] < Epsilon)
					continue;
				Float sx = -normals[i][j].x / normals[i][j].z;
				Float sy = -normals[i][j].y / normals[i][j].z;
				fprintf(fp, "%.8f %.8f\n", sx, sy);
			}
		}
		fclose(fp);
		*/

		return 0;
	}

	Normal sampleNormal(int i, int j, double &weight) {
		ref<Sampler> sampler = m_samplers[Thread::getID() % 233];
		double x = m_aabb.min.x + (j + sampler->next1D()) / (double)m_sqrtSpp * (m_aabb.max.x - m_aabb.min.x);
		double y = m_aabb.min.y + (i + sampler->next1D()) / (double)m_sqrtSpp * (m_aabb.max.y - m_aabb.min.y);
// 		Point o(x, y, 0);
// 		return m_hmap->getNormal(o);

		Point o(x, y, 1e2);
		Ray ray(o, Vector(0, 0, -1.0f), 0);
		Intersection its;
		m_scene->rayIntersect(ray, its);

		Normal normal = its.shFrame.n;

// 		if (normal.z >= 0.9999) {
// 			Log(EInfo, "(%d, %d), h = %.6f", j, i, its.p.z);
// 		}

		int c = math::clamp(math::floorToInt((normal.x + 1.0) * 0.5 * m_size), 0, m_size - 1);
		int r = math::clamp(math::floorToInt((normal.y + 1.0) * 0.5 * m_size), 0, m_size - 1);
		Normal binNormal;
		binNormal.x = (c + 0.5) / (double)m_size * 2.0 - 1.0;
		binNormal.y = (r + 0.5) / (double)m_size * 2.0 - 1.0;
		binNormal.z = std::sqrt(std::max(0.0, 1.0 - binNormal.x * binNormal.x - binNormal.y * binNormal.y));
		if (binNormal.z <= 1e-4)
			Log(EInfo, "%d, %d; (%.6f, %.6f, %.6f)", c, r, binNormal.x, binNormal.y, binNormal.z);

		weight = std::max(0.0, dot(binNormal, m_wi));
		if (weight > 0.0) {
			ray = Ray(its.p + m_wi * ShadowEpsilon, m_wi, 0);
			if (m_maskingOption == 1) {
				m_scene->rayIntersect(ray, its);
				if (its.isValid() && m_aabb.contains(Point2(its.p.x, its.p.y))) {
					weight = 0.0;
					//Log(EInfo, "masked...");
				}
			}
			else if (m_maskingOption == 2) {
				if (m_scene->rayIntersect(ray)) {
					weight = 0.0;
					//Log(EInfo, "masked...");
				}
			}
		}

		return binNormal;
	}

	void outputLeadrNDF(const std::vector<std::vector<Vector> > &normals, char* fname) {
		Vector moments0(0.0);
		Vector moments1(0.0);
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				Vector slope;
				const Vector &norm = normals[i][j];
				if (std::abs(norm.z) < 1e-5)
					continue;
				slope.x = -norm.x / norm.z;
				slope.y = -norm.y / norm.z;

				moments0.x += slope.x;
				moments0.y += slope.y;
				moments1.x += slope.x * slope.x;
				moments1.y += slope.y * slope.y;
				moments1.z += slope.x * slope.y;
			}
		}
		moments0 /= (double)m_spp;
		moments1 /= (double)m_spp;

		Float sigmaX2 = moments1[0] - moments0[0] * moments0[0];
		Float sigmaY2 = moments1[1] - moments0[1] * moments0[1];
		Float cxy = moments1[2] - moments0[0] * moments0[1];
		Float sigmaX = std::sqrt(sigmaX2);
		Float sigmaY = std::sqrt(sigmaY2);
		Float rxy = cxy / (sigmaX * sigmaY);

		Log(EInfo, "LEADR moments computed");
		Log(EInfo, "x = %.6f, y = %.6f", moments0.x, moments0.y);
		Log(EInfo, "x2 = %.6f, y2 = %.6f, xy = %.6f", moments1.x, moments1.y, moments1.z);
		Log(EInfo, "sigma_x = %.6f, sigma_y = %.6f, corr = %.6f", sigmaX, sigmaY, rxy);

		std::vector<std::vector<Vector> > sampledNormals(m_sqrtSpp, std::vector<Vector>(m_sqrtSpp));
#pragma omp parallel for
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				Sampler *sampler = m_samplers[Thread::getID() % 233];
				Point2 sample2 = sampler->next2D();
				Point2 stdSample2 = warp::squareToStdNormal(sample2);		
				Vector2 slope;
				slope.x = stdSample2[0] * sigmaX + moments0[0];
				slope.y = (rxy * stdSample2[0] + std::sqrt(1.0 - rxy * rxy) * stdSample2[1]) * sigmaY + moments0[1];

				Normal wm(-slope.x, -slope.y, 1.0);
				wm = normalize(wm);
				wm /= m_normal_stLim;
				sampledNormals[i][j] = wm;

				//Log(EInfo, "sx = %.6f, sy = %.6f", slope.x, slope.y);
				//Log(EInfo, "%.6f, %.6f, %.6f", wm.x, wm.y, wm.z);
			}
		}

		m_resLEADR = new Bitmap(Bitmap::ELuminance, Bitmap::EFloat32, Vector2i(m_size, m_size));

		Normal avgNormal(0.0);
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				avgNormal += sampledNormals[i][j];
			}
		}
		avgNormal /= (double)m_spp;

		m_resLEADR->clear();
		float *data = m_resLEADR->getFloat32Data();

		double normFactor = (double)(m_size * m_size) / (4.0 * m_spp);
		double dp = 4.0 / (double)(m_size * m_size);
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				const Vector &normal = sampledNormals[i][j];

				//int c = (normal.x > 0.9999f ? m_size - 1 : math::floorToInt((normal.x + 1.0) * 0.5 * m_size));
				//int r = (normal.y > 0.9999f ? m_size - 1 : math::floorToInt((normal.y + 1.0) * 0.5 * m_size));
				int c = math::clamp(math::floorToInt((normal.x + 1.0) * 0.5 * m_size), 0, m_size - 1);
				int r = math::clamp(math::floorToInt((normal.y + 1.0) * 0.5 * m_size), 0, m_size - 1);

				data[(m_size - r - 1) * m_size + c] += normFactor * normal.z;
			}
		}

		// ensure \int D_{wi}(wm)dwm = 1
		double totD = 0.0;
		for (int r = 0; r < m_size; r++) {
			double y = (r + 0.5) / (double)m_size * 2.0 - 1.0;
			for (int c = 0; c < m_size; c++) {
				double x = (c + 0.5) / (double)m_size * 2.0 - 1.0;
				double sinTheta2 = x * x + y * y;
				if (sinTheta2 >= 1.0)
					continue;
				double jacobian = 1.0 / std::sqrt(1.0 - sinTheta2);
				totD += data[r * m_size + c] * jacobian * dp;
			}
		}
		Log(EInfo, "Validation: %.8f should equal to 1", totD);
		
		ref<FileStream> stream = new FileStream(fname, FileStream::ETruncWrite);
		m_resLEADR->write(Bitmap::EOpenEXR, stream);
	}

	void initKMeans(const std::vector<std::vector<Vector> > &normals,
		const std::vector<std::vector<double> > &weights) {
		std::vector<Float> cnt(m_numLobes);
		std::vector<Vector> centers[2];
		
		ref<Sampler> sampler = m_samplers[Thread::getID() % 233];

		// init
		centers[0].resize(m_numLobes);
		centers[1].resize(m_numLobes);
		for (int i = 0; i < m_numLobes; i++) {
			int k;
			while (1) {
				int u = math::floorToInt(sampler->next1D() * m_spp);
				int r = u / m_sqrtSpp;
				int c = u % m_sqrtSpp;
				if (weights[r][c] < Epsilon)
					continue;

				centers[0][i] = normals[r][c];
				break;
			}
		}

		int now = 0;
		int maxIters = 100;
		for (int iter = 0; iter < maxIters; iter++) {
			for (int l = 0; l < m_numLobes; l++) {
				cnt[l] = 0.0;
				centers[1 - now][l] = Vector(0.0);
			}

			for (int i = 0; i < m_sqrtSpp; i++) {
				for (int j = 0; j < m_sqrtSpp; j++) {
					if (weights[i][j] < Epsilon)
						continue;

					// assign cluster
					Float minDist = 1e8;
					int k = -1;
					for (int l = 0; l < m_numLobes; l++) {
						if (centers[now][l].length() < 1e-8)
							continue;
						double d = 1.0 - dot(normals[i][j], normalize(centers[now][l]));
						if (d < minDist) {
							minDist = d;
							k = l;
						}
					}

					// update centers
					Float w = weights[i][j] / normals[i][j].z;
					centers[1 - now][k] += normals[i][j] * w;
					cnt[k] += w;
				}
			}

			now = 1 - now;
			for (int l = 0; l < m_numLobes; l++) {
				if (cnt[l] > 1e-8)
					centers[now][l] /= cnt[l];
			}

			// converged?
			bool converged = true;
			for (int l = 0; l < m_numLobes; l++) {
				if (centers[now][l].length() < 1e-8 || centers[1 - now][l].length() < 1e-8) {
					if (std::abs(centers[now][l].length() - centers[1 - now][l].length()) > 1e-8) {
						converged = false;
						break;
					}
					continue;
				}

				if (dot(normalize(centers[now][l]), normalize(centers[1 - now][l])) < 0.9999) {
					converged = false;
					break;
				}
			}

			if (converged)
				break;
		}

// 		Log(EInfo, "====== kmeans result ======");
// 		for (int l = 0; l < m_numLobes; l++) {
// 			Log(EInfo, "center %d, (%.8f, %.8f, %.8f)", l, centers[now][l].x, centers[now][l].y, centers[now][l].z);
// 		}

		Float totWeights = 0.0;
		m_vmfs = MultiLobeVMF(m_numLobes);
		for (int l = 0; l < m_numLobes; l++)
			totWeights += cnt[l];
		for (int l = 0; l < m_numLobes; l++) {
			m_vmfs.m_alpha[l] = cnt[l] / totWeights;
			if (m_vmfs.m_alpha[l] < 1e-8)
				continue;

			m_vmfs.m_mu[l] = normalize(centers[now][l]);
			Float kappa = VonMisesFisherDistr::forMeanLength(centers[now][l].length());
			kappa = std::min(kappa, 1e4);
			m_vmfs.m_dist[l] = VonMisesFisherDistr(kappa);
		}

// 		Log(EInfo, "====== Init vMF distribution ======");
// 		for (int l = 0; l < m_numLobes; l++) {
// 			Log(EInfo, "alpha = %.6f, kappa = %.6f, mu = (%.6f, %.6f, %.6f)",
// 				m_vmfs.m_alpha[l], m_vmfs.m_dist[l].getKappa(),
// 				m_vmfs.m_mu[l].x, m_vmfs.m_mu[l].y, m_vmfs.m_mu[l].z);
// 		}
	}

	void EM(const std::vector<std::vector<Vector> > &normals,
		const std::vector<std::vector<double> > &weights) {
		std::vector<std::vector<Float> > prob(m_spp, std::vector<Float>(m_numLobes, 0.0));
#pragma omp parallel for
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				if (weights[i][j] < Epsilon)
					continue;
				int idx = i * m_sqrtSpp + j;

				Float sumProb = 0.0;
				for (int l = 0; l < m_numLobes; l++) {
					if (m_vmfs.m_mu[l].length() < 1e-8)
						continue;
					Float cosTheta = dot(m_vmfs.m_mu[l], normals[i][j]);
					prob[idx][l] = m_vmfs.m_dist[l].eval(cosTheta);
					sumProb += prob[idx][l];
				}

				if (sumProb > 1e-8) {
					for (int l = 0; l < m_numLobes; l++)
						prob[idx][l] /= sumProb;
				}
			}
		}

		for (int l = 0; l < m_numLobes; l++) {
			Float alpha = 0.0;
			Vector r(0.0);
			Float totW = 0.0;
			//Float totProb = 0.0;

			for (int i = 0; i < m_sqrtSpp; i++) {
				for (int j = 0; j < m_sqrtSpp; j++) {
					if (weights[i][j] < Epsilon)
						continue;
					int idx = i * m_sqrtSpp + j;

					Float w = weights[i][j] / normals[i][j].z;
					alpha += w * prob[idx][l];
					totW += w;
					
					//totProb += prob[idx][l];
					r += w * prob[idx][l] * normals[i][j];
				}
			}

			if (alpha < 1e-8) {
				m_vmfs.m_mu[l] = Vector(0.0);
				m_vmfs.m_dist[l] = VonMisesFisherDistr(0.0);
				continue;
			}

			r /= alpha;
			//r /= totProb;
			//Log(EInfo, "lobe %d: (%.8f, %.8f, %.8f)", l, r.x, r.y, r.z);
			alpha /= totW;
		
			Float kappa;
			if (std::abs(r.length() - 1.0) < 1e-8) {
				kappa = 1e4;
			}
			else {
				kappa = VonMisesFisherDistr::forMeanLength(r.length());
				kappa = std::min(kappa, 1e4);
			}
			m_vmfs.m_alpha[l] = alpha;
			m_vmfs.m_mu[l] = normalize(r);
			m_vmfs.m_dist[l] = VonMisesFisherDistr(kappa);
		}
	}

	Vector m_wi;
	int m_maskingOption;
	ref<Scene> m_scene;
	Shape *m_hmap;
	ref_vector<Sampler> m_samplers;
	int m_sqrtSpp, m_spp;
	int m_size;
	AABB2 m_aabb;
	double m_normal_stLim;
	
	// D(wm)
	ref<Bitmap> m_D;

	// D_wi(wm)
	ref<Bitmap> m_res;
	
	ref<Bitmap> m_resLEADR;

	int m_numLobes;
	MultiLobeVMF m_vmfs;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(CalcPatchNDF, "Compute the NDF given a patch in the heightmap")
MTS_NAMESPACE_END
