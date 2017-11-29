/*
	Added by Lifan Wu
	Oct 22, 2015
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/render/film.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/util.h>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

class AAFIntegrator : public SamplingIntegrator {
public:
	AAFIntegrator(const Properties &props) : SamplingIntegrator(props) {
		m_initSamples = props.getSize("initSamples", 16);
		m_minSamples = props.getSize("minSamples", 16);
		m_filterRadius = props.getSize("filterRadius", 5);
		m_verbose = props.getBoolean("verbose", false);
		m_adaptive = props.getBoolean("adaptive", true);

		sppMu = 2.f;
		lightSigma = 2.f;

		m_timer = new Timer(false);
	}

	AAFIntegrator(Stream *stream, InstanceManager *manager)
		: SamplingIntegrator(stream, manager) {
		m_subIntegrator = static_cast<SamplingIntegrator *>(manager->getInstance(stream));
		m_initSamples = stream->readInt();
		m_verbose = false;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();

		if (cClass->derivesFrom(MTS_CLASS(Integrator))) {
			if (!cClass->derivesFrom(MTS_CLASS(SamplingIntegrator)))
				Log(EError, "The sub-integrator must be derived from the class SamplingIntegrator");
			m_subIntegrator = static_cast<SamplingIntegrator *>(child);
		}
		else {
			Integrator::addChild(name, child);
		}
	}

	void configureSampler(const Scene *scene, Sampler *sampler) {
		SamplingIntegrator::configureSampler(scene, sampler);
		m_subIntegrator->configureSampler(scene, sampler);
	}

	Float calcOmegaXf(int index) {
		// min{spp_mu / (light_sigma * s2), 1 / (projDist * (1 + s2))}
		return std::min(slopes[0][index].y > 0 ? sppMu / (lightSigma * slopes[0][index].y) : 1e10f, 
			1.f / (projDists[index] * (1.f + slopes[0][index].y)));
	}

	Float calcSpp(Float s1, Float s2, Float wxf, Float projD) {
		//Float spp1 = (1.f / (1.f + s2) + projD * wxf);
		//Float spp2 = (1.f + lightSigma * std::min(s1 * wxf, 1.f / projD * s1 / (1.f + s1)));
		Float spp1 = 1.f + sppMu * s1 / s2;
		Float spp2 = projD / (s2 * lightSigma) + 1.f / (1.f + s2);
		return 4 * spp1 * spp1 * spp2 * spp2;
	}

	Float calcGaussian(Float d, Float sigma) {
		return expf(-d * d / (2.0f * sigma * sigma));
	}

	void filterSlope(Vector2 &slope, bool &occ, int id, int i, int j, int pass) {
		if (i < 0 || i >= width || j < 0 || j >= height)
			return;
		int index = i + j * width;
		if (objId[index] != id)
			return;
		
		if (!occluded[pass][index])
			return;
		slope.x = std::max(slope.x, slopes[pass][index].x);
		slope.y = std::min(slope.y, slopes[pass][index].y);
		occ = (occ || occluded[pass][index]);
	}

	void calcOmegaXf() {
		// filter slopes
		int radius = 5;
		int pixelNum = width * height;

		for (int pass = 0; pass < 2; pass++) {
			for (int x = 0; x < width; x++) {
				for (int y = 0; y < height; y++) {
					int index = x + y * width;
					Vector2 slope = slopes[pass][index];
					bool occ = occluded[pass][index];

					for (int d = -radius; d <= radius; d++) {
						if (pass == 0)
							filterSlope(slope, occ, objId[index], x + d, y, 0);
						else
							filterSlope(slope, occ, objId[index], x, y + d, 1);
					}
					slopes[1 - pass][index] = slope;
					occluded[1 - pass][index] = occ;
				}
			}
		}

		OmegaXf[0].resize(pixelNum);
		OmegaXf[1].resize(pixelNum);
		// calculate beta
		for (int x = 0; x < width; x++) {
			for (int y = 0; y < height; y++) {
				int index = x + y * width;
				OmegaXf[0][index] = calcOmegaXf(index);
			}
		}

		// filter beta
		for (int pass = 0; pass < 2; pass++) {
			for (int x = 0; x < width; x++) {
				for (int y = 0; y < height; y++) {
					int index = x + y * width;
					if (!occluded[0][index])
						continue;

					Float totValue = 0.f;
					Float totWeight = 0.f;
					for (int d = -radius; d <= radius; d++) {
						int i, j;
						if (pass == 0) {
							i = x + d;
							j = y;
						}
						else {
							i = x;
							j = y + d;
						}
						if (i < 0 || i >= width || j < 0 || j >= height)
							continue;
						
						int targetIndex = i + j * width;
						if (!occluded[0][targetIndex] || 
							objId[targetIndex] != objId[index])
							continue;
						Float weight = calcGaussian((Float)d, 0.5f);
						totValue += OmegaXf[pass][targetIndex] * weight;
						totWeight += weight;
					}

					if (totWeight > 1e-4f)
						OmegaXf[1 - pass][index] = totValue / totWeight;
				}
			}
		}

		// calculate spp
		spp.resize(pixelNum);
		for (int x = 0; x < width; x++) {
			for (int y = 0; y < height; y++) {
				int index = x + y * width;
				if (m_adaptive) {
					spp[index] = calcSpp(slopes[0][index].x, slopes[0][index].y, OmegaXf[0][index], projDists[index]);
					spp[index] = math::clamp(spp[index], (Float)m_minSamples, 1024.0);
				}
				else {
					spp[index] = (Float)m_minSamples;
				}
			}
		}
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int sensorResID, int samplerResID) {
		if (!SamplingIntegrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID))
			return false;
		if (m_subIntegrator == NULL)
			Log(EError, "No sub-integrator was specified!");
		Sampler *originalSampler = static_cast<Sampler*>(Scheduler::getInstance()->getResource(samplerResID, 0)); 
		if (originalSampler->getClass()->getName() != "IndependentSampler")
			Log(EError, "The aaf_shadows integrator should only be "
			"used in conjunction with the independent sampler");
	
		Properties props("stratified");
		props.setInteger("sampleCount", m_initSamples);
		props.setInteger("dimension", 8);
		Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
		sampler->configure();

		Sensor *sensor = static_cast<Sensor*>(Scheduler::getInstance()->getResource(sensorResID));

		if (!m_subIntegrator->preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID))
			return false;

		Vector2i filmSize = sensor->getFilm()->getSize();
		bool needsApertureSample = sensor->needsApertureSample();
		bool needsTimeSample = sensor->needsTimeSample();
		const int nSamples = 0;
		Float luminance = 0;

		Point2 apertureSample(0.5f);
		Float timeSample = 0.5f;
		RadianceQueryRecord rRec(scene, sampler);

		int pixelNums = filmSize.x * filmSize.y;
		width = filmSize.x;
		height = filmSize.y;

		positions.resize(pixelNums);
		normals.resize(pixelNums);
		for (int i = 0; i < 2; i++) {
			occluded[i].resize(pixelNums);
			slopes[i].resize(pixelNums);
		}
		objId.resize(pixelNums);
		projDists.resize(pixelNums);
		float *data = new float[pixelNums * 3];

		Float xfov = ((PerspectiveCamera*)sensor)->getXFov();
		m_timer->start();

		for (int x = 0; x < filmSize.x; x++) {
			for (int y = 0; y < filmSize.y; y++) {
				sampler->generate(Point2i(x, y));

				int index = x + filmSize.x * y;
				positions[index] = Point(0.f);
				normals[index] = Vector(0.f);
				Float d2_min = 1e10;
				Float d2_max = 0;
				slopes[0][index].x = slopes[1][index].x = -1e10f;
				slopes[0][index].y = slopes[1][index].y = 1e10f;
				occluded[0][index] = occluded[1][index] = false;
				objId[index] = -1;
				Float cnt = 0.f;

				for (int c = 0; c < m_initSamples; c++) {
					Point2 samplePos(rRec.nextSample2D());
					samplePos.x += x;
					samplePos.y += y;

					if (needsApertureSample)
						apertureSample = rRec.nextSample2D();
					if (needsTimeSample)
						timeSample = rRec.nextSample1D();

					RayDifferential eyeRay;
					Spectrum sampleValue = sensor->sampleRay(
						eyeRay, samplePos, apertureSample, timeSample);

					rRec.newQuery(RadianceQueryRecord::EDistance | RadianceQueryRecord::EIntersection, 
						sensor->getMedium());

					if (rRec.rayIntersect(eyeRay)) {
						positions[index] += rRec.its.p;
						normals[index] += rRec.its.shFrame.n;
						objId[index] = rRec.its.primIndex;
						projDists[index] += 2.f / width * rRec.its.t * tan(xfov * 0.5 * M_PI / 180.f);
						cnt += 1.f;

						Point2 lightSample(rRec.nextSample2D());
						Point2 *sampleArray;
						sampleArray = &lightSample;

						DirectSamplingRecord dRec(rRec.its);

						for (size_t i = 0; i < 1; ++i) {
							scene->sampleEmitterDirect(dRec, sampleArray[i]);
							Float dist2Light;
							//dist2Light = dRec.dist;
							Vector lightDir = dRec.ref - dRec.p;
							lightNormal = dRec.n;
							dist2Light = fabsf(dot(lightDir, dRec.n));
							
							Intersection occ_its_light_to_surface;
							Intersection occ_its_surface_to_light;
							Ray ray_l2s(dRec.p, -dRec.d, Epsilon,
								dRec.dist * (1 - ShadowEpsilon), dRec.time);
							Ray ray_s2l(dRec.ref, dRec.d, Epsilon,
								dRec.dist * (1 - ShadowEpsilon), dRec.time);
							bool hit = scene->rayIntersect(ray_l2s, occ_its_light_to_surface);
							if (hit) {
								occluded[0][index] = true;
								Float t = fabs(dot(occ_its_light_to_surface.t * dRec.d, dRec.n));
								d2_min = std::min(d2_min, t);
								bool hit_s2l = scene->rayIntersect(ray_s2l, occ_its_surface_to_light);
								if (hit_s2l) {
									t = dist2Light - fabs(dot(occ_its_surface_to_light.t * dRec.d, dRec.n));
									d2_max = std::max(d2_max, t);
								}
								else {
									d2_max = std::max(d2_max, dist2Light - t);
								}
								slopes[0][index].x = std::max(slopes[0][index].x, dist2Light / d2_min - 1.f);
								slopes[0][index].y = std::min(slopes[0][index].y, dist2Light / d2_max - 1.f);
							}
						}
					}

					sampler->advance();
				}

				if (cnt > 0) {
					positions[index] /= cnt;
					normals[index] /= cnt;
					if (normals[index].length() > 1e-4)
						normalize(normals[index]);
					projDists[index] /= cnt;
				}
			}
		}

		calcOmegaXf();

		Float preprocessTime = m_timer->getSecondsSinceStart();
		Log(EInfo, "Preprocessing time = %.6fs", preprocessTime);

		for (int y = 0; y < filmSize.y; y++) {
			for (int x = 0; x < filmSize.x; x++) {
				for (int c = 0; c < 3; c++) {
					data[3 * (x + y * filmSize.x) + c] = (occluded[0][x + y * filmSize.x] ? 1.f : 0.f);
				}
			}
		}
		savePfm("occluded.pfm", data, filmSize.x, filmSize.y);

		for (int y = 0; y < filmSize.y; y++) {
			for (int x = 0; x < filmSize.x; x++) {
				int index = x + y * filmSize.x;
				if (!occluded[0][index]) {
					for (int c = 0; c < 3; c++) data[3 * index + c] = 0.f;
					continue;
				}
				data[3 * (x + y * filmSize.x) + 0] = slopes[0][index].x;
				data[3 * (x + y * filmSize.x) + 1] = slopes[0][index].y;
				//data[3 * (x + y * filmSize.x) + 2] = OmegaXf[0][index];
				data[3 * (x + y * filmSize.x) + 2] = spp[index];
			}
		}
		savePfm("s1s2wxf.pfm", data, filmSize.x, filmSize.y);

		return true;
	}

	void renderBlock(const Scene *scene, const Sensor *sensor,
		Sampler *sampler, ImageBlock *block, const bool &stop,
		const std::vector< TPoint2<uint8_t> > &points) const {
		bool needsApertureSample = sensor->needsApertureSample();
		bool needsTimeSample = sensor->needsTimeSample();

		Point2 apertureSample(0.5f);
		Float timeSample = 0.5f;
		RayDifferential sensorRay;

		block->clear();

		uint32_t queryType = RadianceQueryRecord::ESensorRay;

		if (!sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we don't have to */
			queryType &= ~RadianceQueryRecord::EOpacity;

		for (size_t i = 0; i < points.size(); ++i) {
			Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
			if (stop)
				break;

			int index = offset.x + offset.y * width;
			int sqrtSpp = (int)(ceil(sqrtf(spp[index])));

			RadianceQueryRecord rRec(scene, sampler);
			Float diffScaleFactor = 1.0f / (Float)sqrtSpp;

			sampler->generate(offset);

			for (size_t j = 0; j < sqrtSpp * sqrtSpp; j++) {
				rRec.newQuery(queryType, sensor->getMedium());

				Point2 samplePos(Vector2(rRec.nextSample2D()));
				int x = j % sqrtSpp;
				int y = j / sqrtSpp;
				samplePos.x = ((Float)x + samplePos.x) / (Float)sqrtSpp + offset.x;
				samplePos.y = ((Float)y + samplePos.y) / (Float)sqrtSpp + offset.y;

				if (needsApertureSample)
					apertureSample = rRec.nextSample2D();
				if (needsTimeSample)
					timeSample = rRec.nextSample1D();

				Spectrum spec = sensor->sampleRayDifferential(
					sensorRay, samplePos, apertureSample, timeSample);

				sensorRay.scaleDifferential(diffScaleFactor);

				spec *= m_subIntegrator->Li(sensorRay, rRec);
				block->put(samplePos, spec, rRec.alpha);
				
				sampler->advance();
			}
		}
	}

	Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
		return m_subIntegrator->Li(ray, rRec);
	}

	Spectrum E(const Scene *scene, const Intersection &its, const Medium *medium,
		Sampler *sampler, int nSamples, bool includeIndirect) const {
		return m_subIntegrator->E(scene, its, medium,
			sampler, nSamples, includeIndirect);
	}

	void postprocess(const Scene *scene, RenderQueue *queue,
		const RenderJob *job, int sceneResID, int sensorResID,
		int samplerResID) {
		const Film *film = scene->getSensor()->getFilm();
		const Bitmap *bitmap = film->getImageBlock()->getBitmap();
		ref<Bitmap> img = new Bitmap(*bitmap);

		int pixelNum = height * width;
		imgRes[0].resize(pixelNum);
		imgRes[1].resize(pixelNum);

		//totWeights[0].resize(pixelNum);
		//totWeights[1].resize(pixelNum);

		for (int x = 0; x < width; x++) {
			for (int y = 0; y < height; y++) {
				int index = x + y * width;
				imgRes[0][index] = img->getPixel(Point2i(x, y));
			}
		}

		float *data = new float[pixelNum * 3];
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int index = x + y * width;
				Spectrum color(imgRes[0][index]);
				for (int c = 0; c < 3; c++) {
					data[3 * (x + y * width) + c] = color[c];
				}
			}
		}
		savePfm("origin.pfm", data, width, height);

		const float dist_scale_threshold = 0.5f;
		const float dist_threshold = 1.f;
		const float cos_angle_threshold = 0.7f;

		for (int pass = 0; pass < 2; pass++) {
			for (int x = 0; x < width; x++) {
				for (int y = 0; y < height; y++) {
					int index = x + y * width;
					if (!occluded[0][index])
						continue;

					Spectrum totValue(0.f);
					Float totWeight = 0.f;
					for (int d = -m_filterRadius; d <= m_filterRadius; d++) {
						int i, j;
						if (pass == 0) {
							i = x + d;
							j = y;
						}
						else {
							i = x;
							j = y + d;
						}
						if (i < 0 || i >= width || j < 0 || j >= height)
							continue;

						int targetIndex = i + j * width;
						if (!occluded[0][targetIndex] ||
							objId[targetIndex] != objId[index])
							continue;
						if (fabsf(OmegaXf[0][targetIndex] - OmegaXf[0][index]) > dist_scale_threshold)
							continue;
						
						Vector diff = positions[targetIndex] - positions[index];
						Float normComp = dot(diff, lightNormal);
						Float dist2 = diff.lengthSquared() - normComp * normComp;
						if (fabsf(dist2) > dist_threshold)
							continue;

						if (fabsf(dot(normals[targetIndex], normals[index])) <= cos_angle_threshold)
							continue;

						Float wxf = OmegaXf[0][index];
						Float weight = expf(-9.f * dist2 * wxf * wxf * 0.5f);
						totValue += imgRes[pass][targetIndex] * weight;
						totWeight += weight;
					}

					if (totWeight > 1e-4f) {
						imgRes[1 - pass][index] = totValue / totWeight;
						//totWeights[pass][index] = totWeight;
					}
				}
			}
		}

// 		for (int y = 0; y < height; y++) {
// 			for (int x = 0; x < width; x++) {
// 				int index = x + y * width;
// 				data[3 * (x + y * width) + 0] = totWeights[0][index];
// 				data[3 * (x + y * width) + 1] = totWeights[1][index];
// 				data[3 * (x + y * width) + 2] = 0;
// 			}
// 		}
// 		savePfm("totWeights.pfm", data, width, height);

		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int index = x + y * width;
				Spectrum color(imgRes[0][index]);
				for (int c = 0; c < 3; c++) {
					data[3 * (x + y * width) + c] = color[c];
				}
			}
		}
		savePfm("filtered.pfm", data, width, height);

		/*
		img = img->convert(Bitmap::ERGB, Bitmap::EFloat32);

		fs::path output = scene->getDestinationFile();
		output.replace_extension("test.pfm");
		Log(EInfo, "output path = %s", output.string().c_str());
		ref<FileStream> stream = new FileStream(output, FileStream::ETruncWrite);

		img->write(Bitmap::EPFM, stream);
		*/
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		SamplingIntegrator::serialize(stream, manager);
		manager->serialize(stream, m_subIntegrator.get());
		
		stream->writeInt(m_initSamples);
	}

	void bindUsedResources(ParallelProcess *proc) const {
		m_subIntegrator->bindUsedResources(proc);
	}

	void wakeup(ConfigurableObject *parent,
		std::map<std::string, SerializableObject *> &params) {
		m_subIntegrator->wakeup(this, params);
	}

	void cancel() {
		SamplingIntegrator::cancel();
		m_subIntegrator->cancel();
	}

	const Integrator *getSubIntegrator(int idx) const {
		if (idx != 0)
			return NULL;
		return m_subIntegrator.get();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "AdaptiveIntegrator[" << endl
			<< "  initSamples = " << m_initSamples << "," << endl
			<< "  subIntegrator = " << indent(m_subIntegrator->toString()) << endl
			<< "]";
		return oss.str();
	}

	void savePfm(char *fileName, float *data, int width, int height) {
		char header[512];
		sprintf(header, "PF\n%d %d\n-1.000000\n", width, height);

		FILE *fp = fopen(fileName, "wb");
		fwrite(header, strlen(header), 1, fp);
		for (int y = height - 1; y >= 0; y--) {
			for (int x = 0; x < width; x++) {
				fwrite(&data[3 * (y * width + x) + 0], sizeof(float), 1, fp);
				fwrite(&data[3 * (y * width + x) + 1], sizeof(float), 1, fp);
				fwrite(&data[3 * (y * width + x) + 2], sizeof(float), 1, fp);
			}
		}
		fclose(fp);
	}

	MTS_DECLARE_CLASS()
private:
	ref<SamplingIntegrator> m_subIntegrator;
	int m_initSamples;
	int m_minSamples;
	int m_filterRadius;
	bool m_adaptive;

	Timer *m_timer;

	int width, height;
	
	// 0: max_slope, 1: min_slope
	std::vector<Vector2> slopes[2];
	std::vector<bool> occluded[2];
	
	std::vector<Vector> normals;
	std::vector<Point> positions;
	std::vector<int> objId;
	std::vector<Float> projDists;

	std::vector<Float> OmegaXf[2];
	std::vector<Float> spp;
	std::vector<Spectrum> imgRes[2];

	//std::vector<Float> totWeights[2];

	// assume only one area light source
	Vector lightNormal;
	Float lightSigma;
	Float sppMu;

	bool m_verbose;
};

MTS_IMPLEMENT_CLASS_S(AAFIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(AAFIntegrator, "Axis-aligned filtering based integrator");
MTS_NAMESPACE_END
