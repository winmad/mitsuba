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

		double normFactor = (double)(m_size * m_size) / (4.0 * m_spp);
		
#pragma omp parallel for
		for (int i = 0; i < m_sqrtSpp; i++) {
			for (int j = 0; j < m_sqrtSpp; j++) {
				normals[i][j] = sampleNormal(i, j, weights[i][j]);
				normals[i][j] /= m_normal_stLim;
			}
		}
		Log(EInfo, "Finish binning.");

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
				totD += data[r * m_size + c] * jacobian * dp;
			}
		}
		Log(EInfo, "Validation: %.8f should equal to 1", totD);

		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		m_res->write(Bitmap::EOpenEXR, stream);

		// fitting with multi-lobe vMF
		if (argc > 14) {
			m_numLobes = std::atoi(argv[14]);
			initKMeans(normals, weights);

			char fname[256];
			sprintf(fname, "vmf_init_%s", argv[13]);
			m_vmfs.outputDistribution(m_size, fname);

			int maxIters = 10;
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
		}

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

		weight = std::max(0.0, dot(normal, m_wi));
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

		return normal;
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
	ref<Bitmap> m_res;

	int m_numLobes;
	MultiLobeVMF m_vmfs;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(CalcPatchNDF, "Compute the NDF given a patch in the heightmap")
MTS_NAMESPACE_END