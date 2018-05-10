#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/core/aabb.h>
#include <mitsuba/render/range.h>
#include <mitsuba/render/spherical_distribution.h>
#include <boost/filesystem/path.hpp>

MTS_NAMESPACE_BEGIN

class EffBRDFDirect : public Utility {
public:
	int run(int argc, char **argv) {
		m_scene = loadScene(argv[1]);
		m_xmin = std::atof(argv[2]);
		m_xmax = std::atof(argv[3]);
		m_ymin = std::atof(argv[4]);
		m_ymax = std::atof(argv[5]);
		m_sqrtSpp = std::atoi(argv[6]);

		m_wiResolution = std::atoi(argv[7]);
		m_woResolution = std::atoi(argv[8]);
		m_shadowOption = std::atoi(argv[9]);
		char *filename = argv[10];

		m_aabb = AABB2(Point2(m_xmin, m_ymin), Point2(m_xmax, m_ymax));
		m_scene->initialize();
		m_spp = m_sqrtSpp * m_sqrtSpp;

		// init samplers
		Properties props = Properties("independent");
		props.setInteger("seed", 19931004);
		m_samplers.resize(233);
		m_samplers[0] = static_cast<Sampler *>(PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));
		m_samplers[0]->configure();
		for (int i = 1; i < 233; i++) {
			m_samplers[i] = m_samplers[0]->clone();
		}

		int numRows = m_wiResolution * m_wiResolution;
		int numCols = m_woResolution * m_woResolution;
		m_effBRDF.resize(numRows);
		for (int i = 0; i < numRows; i++) {
			m_effBRDF[i].resize(numCols);
			for (int j = 0; j < numCols; j++) {
				m_effBRDF[i][j] = Spectrum(0.0);
			}
		}

		Float smallest = 0.02;
		genLinspace(0.0, 1.0, m_wiResolution, smallest, m_wiVec);
		genLinspace(0.0, 1.0, m_woResolution, smallest, m_woVec);

#pragma omp parallel for
		for (int wiIdx = 0; wiIdx < m_wiResolution * m_wiResolution; wiIdx++) {
			int r1 = wiIdx / m_wiResolution;
			int c1 = wiIdx % m_wiResolution;
			Vector wi = warp::squareToUniformHemisphereConcentric(Point2(m_wiVec[c1], m_wiVec[r1]));
			Log(EInfo, "filling wiIdx = %d", wiIdx);

			for (int woIdx = 0; woIdx < m_woResolution * m_woResolution; woIdx++) {
				int r2 = woIdx / m_woResolution;
				int c2 = woIdx % m_woResolution;
				Vector wo = warp::squareToUniformHemisphereConcentric(Point2(m_woVec[c2], m_woVec[r2]));

				m_effBRDF[wiIdx][woIdx] = stocEvalBRDF(wi, wo, m_samplers[Thread::getID() % 233]);
			}
		}

		Log(EInfo, "Finished. Output 4D effective BRDF");
		output4DEffBRDF(filename);

// 		int tr = 12;
// 		int tc = 12;
// 		Vector wi = warp::squareToUniformHemisphereConcentric(Point2(m_wiVec[tc], m_wiVec[tr]));
// 		//Vector wi(0.9969, 0.0, 0.0784);
// 		//Vector wi(0.987162, 0.0, 0.1597222);
// 		//Vector wi(0.95217425, 0.0, 0.3055556);
// 		//wi = normalize(wi);
// 		output2DSlice(wi, filename);
	}

	void genLinspace(Float st, Float ed, int n, Float smallest, std::vector<Float> &vec) {
		vec.resize(n);
		Float step = (ed - st) / (double)(n - 1);
		for (int i = 0; i < n; i++) {
			vec[i] = st + step * i;
		}
		vec[0] = smallest;
		vec[n - 1] = ed - smallest;
	}

	Spectrum stocEvalBRDF(const Vector &wi, const Vector &wo, Sampler *sampler) {
		Spectrum res(0.0);
		Float totWeight = 0.0;
		
		Intersection its;
		Intersection shadowIts;
		for (int r = 0; r < m_sqrtSpp; r++) {
			double y = m_aabb.min.y + (r + sampler->next1D()) / (double)m_sqrtSpp * (m_aabb.max.y - m_aabb.min.y);
			for (int c = 0; c < m_sqrtSpp; c++) {
				double x = m_aabb.min.x + (c + sampler->next1D()) / (double)m_sqrtSpp * (m_aabb.max.x - m_aabb.min.x);

				// get point on the micro-surface
				Point o(x, y, 1e2);
				Ray ray(o, Vector(0, 0, -1.0), 0);
				
				m_scene->rayIntersect(ray, its);
				
				Normal normal = its.shFrame.n;
				double cosG = std::max(0.0, normal.z);
				if (cosG < 1e-5)
					continue;

				// masking: visible projected area
				Float weight = std::max(0.0, dot(normal, wi)) / cosG;
				if (weight > 0.0) {
					ray = Ray(its.p + wi * ShadowEpsilon, wi, 0);
					if (m_shadowOption == 1) {
						m_scene->rayIntersect(ray, shadowIts);
						if (shadowIts.isValid() && m_aabb.contains(Point2(shadowIts.p.x, shadowIts.p.y))) {
							weight = 0.0;
						}
					} else if (m_shadowOption == 2) {
						if (m_scene->rayIntersect(ray)) {
							weight = 0.0;
						}
					}
				}
				if (weight < 1e-5)
					continue;
				totWeight += weight;

				// shadowing
				Float cosWo = std::max(0.0, dot(normal, wo));
				if (cosWo < 1e-5)
					continue;
				ray = Ray(its.p + wo * ShadowEpsilon, wo, 0);
				if (m_shadowOption == 1) {
					m_scene->rayIntersect(ray, shadowIts);
					if (shadowIts.isValid() && m_aabb.contains(Point2(shadowIts.p.x, shadowIts.p.y))) {
						cosWo = 0.0;
					}
				} else if (m_shadowOption == 2) {
					if (m_scene->rayIntersect(ray)) {
						cosWo = 0.0;
					}
				}
				if (cosWo < 1e-5)
					continue;

				// eval
				BSDFSamplingRecord bRec(its, its.toLocal(wi), its.toLocal(wo));
				Spectrum spec = its.getBSDF()->eval(bRec);
				res += spec * cosWo * weight;
			}
		}

		if (totWeight > 1e-5)
			res /= totWeight;
		return res;
	}

	void output2DSlice(Vector wi, char *filename) {
		ref<Bitmap> lobe = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, Vector2i(m_woResolution));
		float *data = lobe->getFloat32Data();

#pragma omp parallel for
		for (int r2 = 0; r2 < m_woResolution; r2++) {
			for (int c2 = 0; c2 < m_woResolution; c2++) {
				Vector wo = warp::squareToUniformHemisphereConcentric(Point2(m_woVec[c2], m_woVec[r2]));
				Spectrum res = stocEvalBRDF(wi, wo, m_samplers[Thread::getID() % 233]);

				int woIdx = r2 * m_woResolution + c2;
				for (int k = 0; k < 3; k++)
					data[3 * woIdx + k] = res[k];
			}
		}
		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		lobe->write(Bitmap::EOpenEXR, stream);
	}

	void output4DEffBRDF(char *filename) const {
		ref<Bitmap> bitmap = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, 
			Vector2i(m_woResolution * m_woResolution, m_wiResolution * m_wiResolution));
		float *data = bitmap->getFloat32Data();
		for (int r = 0; r < m_wiResolution * m_wiResolution; r++) {
			for (int c = 0; c < m_woResolution * m_woResolution; c++) {
				*data++ = m_effBRDF[r][c][0];
				*data++ = m_effBRDF[r][c][1];
				*data++ = m_effBRDF[r][c][2];
			}
		}
		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		bitmap->write(Bitmap::EOpenEXR, stream);
	}

	ref<Scene> m_scene;
	ref_vector<Sampler> m_samplers;
	int m_sqrtSpp, m_spp;
	int m_wiResolution, m_woResolution;
	double m_xmin, m_xmax, m_ymin, m_ymax;
	AABB2 m_aabb;

	// 1: local, 2: global
	int m_shadowOption;

	std::vector<std::vector<Spectrum> > m_effBRDF;
	std::vector<Float> m_wiVec, m_woVec;

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(EffBRDFDirect, "Evaluate effective BRDF (direct illumination)")
MTS_NAMESPACE_END