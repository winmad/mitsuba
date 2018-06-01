#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <boost/filesystem/path.hpp>
#include "multi_vmf.h"

MTS_NAMESPACE_BEGIN

class NormalMipmap : public Utility {
public:
	typedef std::vector<std::vector<Vector> > Vector2DArray;

	int run(int argc, char **argv) {
		m_numLobes = std::atoi(argv[1]);
		m_resolution = std::atoi(argv[2]);
		m_sqrtSpp = std::atoi(argv[3]);
		m_maxScale = std::atoi(argv[4]);
		m_scene = loadScene(argv[5]);

		m_scene->initialize();
		m_spp = m_sqrtSpp * m_sqrtSpp;
		m_hmap = m_scene->getShapes()[0];

		m_samplers.resize(233);
		m_samplers[0] = static_cast<Sampler *>(PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));
		m_samplers[0]->configure();
		for (int i = 1; i < 233; i++) {
			m_samplers[i] = m_samplers[0]->clone();
		}

		// sample normals
		m_normals = std::vector<std::vector<Vector> >(m_resolution * m_sqrtSpp,
			std::vector<Vector>(m_resolution * m_sqrtSpp));
//#pragma omp parallel for
		for (int i = 0; i < m_resolution * m_sqrtSpp; i++) {
			for (int j = 0; j < m_resolution * m_sqrtSpp; j++) {
				Point2 uv;
				uv.x = (Float)(j + m_samplers[0]->next1D()) / (Float)(m_resolution * m_sqrtSpp);
				uv.y = (Float)(i + m_samplers[0]->next1D()) / (Float)(m_resolution * m_sqrtSpp);

				Normal norm;
				m_hmap->getPosAndNormal(uv, NULL, &norm);
				m_normals[i][j] = norm;
			}
		}
		Log(EInfo, "Finish sampling normals...");

		//for (int i = 2; i <= m_maxScale; i *= 2)
		//	generateVMFs(i);
		generateVMFs(m_maxScale);

		return 0;
	}

	void generateVMFs(int scale) {
		std::vector<ref<Bitmap> > bitmaps(m_numLobes);
		for (int l = 0; l < m_numLobes; l++) {
			bitmaps[l] = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, 
				Vector2i(m_resolution / scale, m_resolution / scale * 2));
		}

		int step = scale * m_sqrtSpp;
		char fname[256];

#pragma omp parallel for
		for (int i = 0; i < m_resolution * m_sqrtSpp; i += step) {
			for (int j = 0; j < m_resolution * m_sqrtSpp; j += step) {
				MultiLobeVMF vmfs;
				initKMeans(m_normals, i, step, j, step, vmfs);
				
// 				sprintf(fname, "vmf_init_%d_%d.exr", i / step, j / step);
// 				vmfs.outputDistribution(128, fname);

				int maxIters = 20;
				for (int iter = 0; iter < maxIters; iter++) {
					EM(m_normals, i, step, j, step, vmfs);
				}

// 				sprintf(fname, "vmf_%d_%d.exr", i / step, j / step);
// 				vmfs.outputDistribution(128, fname);
// 
// 				sprintf(fname, "ndf_%d_%d.exr", i / step, j / step);
// 				outputOriginalNDF(m_normals, i, step, j, step, fname);

				putVMF(vmfs, bitmaps, j / step, i / step, m_resolution / scale,
					m_resolution / scale);

				// no lobe error
				bool flag = false;
				for (int l = 0; l < m_numLobes; l++) {
					if (vmfs.m_alpha[l] > 1e-5)
						flag = true;
				}
				if (!flag) {
					std::ostringstream oss;

					/*
					for (int di = 0; di < step; di++) {
						for (int dj = 0; dj < step; dj++) {
							Vector &n = m_normals[i + di][j + dj];
							oss << "(" << n.x << ", " << n.y << ", " << n.z << ")\n";
						}
					}
					*/

					initKMeans(m_normals, i, step, j, step, vmfs);
					for (int l = 0; l < m_numLobes; l++) {
						oss << "alpha: " << vmfs.m_alpha[l] << " kappa: " << vmfs.m_dist[l].getKappa() << "\n";
					}

					Log(EInfo, "%s", oss.str().c_str());					
					exit(0);
				}

// 				if (i / step == 0 && j / step == 0) {
// 					sprintf(fname, "vmf_%d_%d.exr", i / step, j / step);
// 					vmfs.outputDistribution(128, fname);
// 					 
// 					sprintf(fname, "ndf_%d_%d.exr", i / step, j / step);
// 					outputOriginalNDF(m_normals, i, step, j, step, fname);
// 				}
			}
		}

		for (int l = 0; l < m_numLobes; l++) {
			char fname[256];
			sprintf(fname, "vmf_%dx_lobe_%d.exr", scale, l);
			ref<FileStream> stream = new FileStream(fs::path(fname), FileStream::ETruncWrite);
			bitmaps[l]->write(Bitmap::EOpenEXR, stream);
		}
	}

	void initKMeans(const Vector2DArray &normals, int rSt, int rSize, int cSt, int cSize,
		MultiLobeVMF &vmfsInit) {
		std::vector<Float> cnt(m_numLobes);
		std::vector<Vector> centers[2];

		ref<Sampler> sampler = m_samplers[Thread::getID() % 233];

		// init
		centers[0].resize(m_numLobes);
		centers[1].resize(m_numLobes);
		for (int i = 0; i < m_numLobes; i++) {
			int counter = 0;
			while (1) {
				int u = math::floorToInt(sampler->next1D() * rSize);
				int v = math::floorToInt(sampler->next1D() * cSize);
				int r = rSt + v;
				int c = cSt + u;

				if (counter < 1000) {
					counter++;
					for (int l = 0; l < i; l++) {
						if (dot(normals[r][c], centers[0][l]) > 0.9999)
							continue;
					}
				}
			
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

			for (int i = rSt; i < rSt + rSize; i++) {
				for (int j = cSt; j < cSt + cSize; j++) {
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
					centers[1 - now][k] += normals[i][j];
					cnt[k] += 1.0;
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
		vmfsInit = MultiLobeVMF(m_numLobes);
		for (int l = 0; l < m_numLobes; l++)
			totWeights += cnt[l];
		for (int l = 0; l < m_numLobes; l++) {
			vmfsInit.m_alpha[l] = cnt[l] / totWeights;
			if (vmfsInit.m_alpha[l] < 1e-8)
				continue;

			vmfsInit.m_mu[l] = normalize(centers[now][l]);

			Float kappa;
			if (std::abs(centers[now][l].length() - 1.0) < 1e-3)
				kappa = 1e3;
			else {
				kappa = VonMisesFisherDistr::forMeanLength(centers[now][l].length());
				kappa = std::min(kappa, 1e3);
			}
			vmfsInit.m_dist[l] = VonMisesFisherDistr(kappa);
		}
	}

	void EM(const Vector2DArray &normals, int rSt, int rSize, int cSt, int cSize,
		MultiLobeVMF &vmfs) {
		std::vector<std::vector<Float> > prob(rSize * cSize, std::vector<Float>(m_numLobes, 0.0));
#pragma omp parallel for
		for (int i = 0; i < rSize; i++) {
			for (int j = 0; j < cSize; j++) {
				int idx = i * cSize + j;

				Float sumProb = 0.0;
				for (int l = 0; l < m_numLobes; l++) {
					if (vmfs.m_mu[l].length() < 1e-8)
						continue;
					Float cosTheta = dot(vmfs.m_mu[l], normals[rSt + i][cSt + j]);
					prob[idx][l] = vmfs.m_dist[l].eval(cosTheta);
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

			for (int i = 0; i < rSize; i++) {
				for (int j = 0; j < cSize; j++) {
					int idx = i * cSize + j;

					Float w = 1.0;
					alpha += w * prob[idx][l];
					totW += w;

					r += w * prob[idx][l] * normals[rSt + i][cSt + j];
				}
			}

			if (alpha < 1e-8) {
				vmfs.m_alpha[l] = 0.0;
				vmfs.m_mu[l] = Vector(0.0);
				vmfs.m_dist[l] = VonMisesFisherDistr(0.0);
				continue;
			}

			r /= alpha;
			alpha /= totW;

			Float kappa;
			if (std::abs(r.length() - 1.0) < 1e-8) {
				kappa = 1e3;
			}
			else {
				kappa = VonMisesFisherDistr::forMeanLength(r.length());
				kappa = std::min(kappa, 1e3);
			}
			vmfs.m_alpha[l] = alpha;
			vmfs.m_mu[l] = normalize(r);
			vmfs.m_dist[l] = VonMisesFisherDistr(kappa);
		}
	}

	void putVMF(const MultiLobeVMF &vmfs, std::vector<ref<Bitmap> > &bitmaps, int x, int y, 
		int width, int heightOffset) {
		for (int l = 0; l < m_numLobes; l++) {
			float *data = bitmaps[l]->getFloat32Data();
			int idx = x + y * width;
			data[3 * idx + 0] = vmfs.m_alpha[l];
			data[3 * idx + 1] = vmfs.m_dist[l].getKappa();
			data[3 * idx + 2] = 0.f;

			idx = x + (y + heightOffset) * width;
			data[3 * idx + 0] = 0.5 * (vmfs.m_mu[l].x + 1.0);
			data[3 * idx + 1] = 0.5 * (vmfs.m_mu[l].y + 1.0);
			data[3 * idx + 2] = 0.5 * (vmfs.m_mu[l].z + 1.0);

			if (vmfs.m_alpha[l] > 1e-8 && vmfs.m_mu[l].length() < 1e-8) {
				Log(EInfo, "%.8f, %.8f, (%.8f, %.8f, %.8f)",
					vmfs.m_alpha[l], vmfs.m_dist[l].getKappa(),
					vmfs.m_mu[l].x, vmfs.m_mu[l].y, vmfs.m_mu[l].z);
			}
		}
	}

	void outputOriginalNDF(const Vector2DArray &normals, int rSt, int rSize, int cSt, int cSize,
		char *filename) {
		int size = 128;
		ref<Bitmap> res = new Bitmap(Bitmap::ELuminance, Bitmap::EFloat32, Vector2i(size, size));
		res->clear();
		float *data = res->getFloat32Data();

		for (int i = rSt; i < rSt + rSize; i++) {
			for (int j = cSt; j < cSt + cSize; j++) {
				const Vector &normal = normals[i][j];

				Point2 p = warp::uniformHemisphereToSquareConcentric(normal);
				int c = math::clamp(math::floorToInt(p.x * size), 0, size - 1);
				int r = math::clamp(math::floorToInt(p.y * size), 0, size - 1);

				data[r * size + c] += size * size / (2.0 * M_PI * rSize * cSize);
			}
		}

		ref<FileStream> stream = new FileStream(fs::path(filename), FileStream::ETruncWrite);
		res->write(Bitmap::EOpenEXR, stream);
	}

	int m_numLobes;
	int m_resolution;
	int m_sqrtSpp, m_spp;
	int m_maxScale;
	ref<Scene> m_scene;
	Shape *m_hmap;
	ref_vector<Sampler> m_samplers;

	Vector2DArray m_normals;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(NormalMipmap, "Generate top-visible normal mip-mapping")
MTS_NAMESPACE_END
