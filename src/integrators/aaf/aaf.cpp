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
		m_filterRadius = props.getSize("filterRadius", 5);
		m_verbose = props.getBoolean("verbose", false);

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

	void calcFilterBeta() {
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

		// calculate beta


		// filter beta
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int sensorResID, int samplerResID) {
		if (!SamplingIntegrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID))
			return false;
		if (m_subIntegrator == NULL)
			Log(EError, "No sub-integrator was specified!");
//		Sampler *sampler = static_cast<Sampler *>(Scheduler::getInstance()->getResource(samplerResID, 0)); 
		Properties props("stratified");
		props.setInteger("sampleCount", m_initSamples);
		props.setInteger("dimension", 7);
		Sampler *sampler = static_cast<Sampler*> (PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
		sampler->configure();

		Sensor *sensor = static_cast<Sensor *>(Scheduler::getInstance()->getResource(sensorResID));
// 		if (sampler->getClass()->getName() != "IndependentSampler")
// 			Log(EError, "The error-controlling integrator should only be "
// 			"used in conjunction with the independent sampler");
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
		dists.resize(pixelNums);
		for (int i = 0; i < 2; i++) {
			occluded[i].resize(pixelNums);
			slopes[i].resize(pixelNums);
		}
		objId.resize(pixelNums);
		Float *data = new Float[pixelNums * 3];

		m_timer->start();

		for (int x = 0; x < filmSize.x; x++) {
			for (int y = 0; y < filmSize.y; y++) {
				sampler->generate(Point2i(x, y));

				int index = x + filmSize.x * y;
				positions[index] = Point(0.f);
				normals[index] = Vector(0.f);
				dists[index].x = dists[index].z = 1e10f;
				dists[index].y = dists[index].w = 0;
				slopes[0][index].x = slopes[1][index].x = -1e10f;
				slopes[0][index].y = slopes[1][index].y = 1e10f;
				occluded[0][index] = occluded[1][index] = false;
				objId[index] = -1;

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

						Point2 lightSample(rRec.nextSample2D());
						Point2 *sampleArray;
						sampleArray = &lightSample;

						DirectSamplingRecord dRec(rRec.its);

						for (size_t i = 0; i < 1; ++i) {
							scene->sampleEmitterDirect(dRec, sampleArray[i]);
							Float dist2Light = dRec.dist;
							
							Intersection occ_its_light_to_surface;
							Intersection occ_its_surface_to_light;
							Ray ray_l2s(dRec.p, -dRec.d, Epsilon,
								dRec.dist * (1 - ShadowEpsilon), dRec.time);
							Ray ray_s2l(dRec.ref, dRec.d, Epsilon,
								dRec.dist * (1 - ShadowEpsilon), dRec.time);
							bool hit = scene->rayIntersect(ray_l2s, occ_its_light_to_surface);
							if (hit) {
								occluded[index] = true;
								dists[index].x = std::min(dists[index].x, dist2Light);
								dists[index].y = std::max(dists[index].y, dist2Light);
								dists[index].z = std::min(dists[index].z, occ_its_light_to_surface.t);
								bool hit_s2l = scene->rayIntersect(ray_s2l, occ_its_surface_to_light);
								Assert(hit_s2l);
								dists[index].w = std::max(dists[index].w, dist2Light - occ_its_surface_to_light.t);
							}
						}
					}

					if (occluded[index]) {
						slopes[index].x = dists[index].y / dists[index].z - 1.f;
						slopes[index].y = dists[index].x / dists[index].w - 1.f;
					}

					sampler->advance();
				}

				positions[index] /= (Float)m_initSamples;
				normals[index] /= (Float)m_initSamples;
				if (normals[index].length() > 1e-4)
					normalize(normals[index]);
			}
		}

		Float preprocessTime = m_timer->getSecondsSinceStart();
		Log(EInfo, "Preprocessing time = %.6fs", preprocessTime);

		for (int y = 0; y < filmSize.y; y++) {
			for (int x = 0; x < filmSize.x; x++) {
				for (int c = 0; c < 3; c++) {
					data[3 * (x + y * filmSize.x) + c] = (occluded[x + y * filmSize.x] ? 1.f : 0.f);
				}
			}
		}
		savePfm("occluded.pfm", data, filmSize.x, filmSize.y);

		for (int y = 0; y < filmSize.y; y++) {
			for (int x = 0; x < filmSize.x; x++) {
				int index = x + y * filmSize.x;
				if (!occluded[index]) {
					for (int c = 0; c < 3; c++) data[3 * index + c] = 0.f;
					continue;
				}
				//Spectrum s1 = heatMap(slopes[index].x);
				Spectrum s1(slopes[index].x);
				for (int c = 0; c < 3; c++) {
					data[3 * (x + y * filmSize.x) + c] = s1[c];
				}
			}
		}
		savePfm("s1.pfm", data, filmSize.x, filmSize.y);

		for (int y = 0; y < filmSize.y; y++) {
			for (int x = 0; x < filmSize.x; x++) {
				int index = x + y * filmSize.x;
				if (!occluded[index]) {
					for (int c = 0; c < 3; c++) data[3 * index + c] = 0.f;
					continue;
				}
				//Spectrum s2 = heatMap(slopes[index].y);
				Spectrum s2(slopes[index].y);
				for (int c = 0; c < 3; c++) {
					data[3 * (x + y * filmSize.x) + c] = s2[c];
				}
			}
		}
		savePfm("s2.pfm", data, filmSize.x, filmSize.y);

		return true;
	}

	void renderBlock(const Scene *scene, const Sensor *sensor,
		Sampler *sampler, ImageBlock *block, const bool &stop,
		const std::vector< TPoint2<uint8_t> > &points) const {
		Float diffScaleFactor = 1.0f /
			std::sqrt((Float)sampler->getSampleCount());

		bool needsApertureSample = sensor->needsApertureSample();
		bool needsTimeSample = sensor->needsTimeSample();

		RadianceQueryRecord rRec(scene, sampler);
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

			sampler->generate(offset);

			for (size_t j = 0; j < sampler->getSampleCount(); j++) {
				rRec.newQuery(queryType, sensor->getMedium());
				Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

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
		img = img->convert(Bitmap::ERGB, Bitmap::EFloat32);

		fs::path output = scene->getDestinationFile();
		output.replace_extension("test.pfm");
		Log(EInfo, "output path = %s", output.string().c_str());
		ref<FileStream> stream = new FileStream(output, FileStream::ETruncWrite);

		img->write(Bitmap::EPFM, stream);
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
				fwrite(&data[3 * (y * width + x) + 2], sizeof(float), 1, fp);
				fwrite(&data[3 * (y * width + x) + 1], sizeof(float), 1, fp);
				fwrite(&data[3 * (y * width + x) + 0], sizeof(float), 1, fp);
			}
		}
		fclose(fp);
	}

	MTS_DECLARE_CLASS()
private:
	ref<SamplingIntegrator> m_subIntegrator;
	int m_initSamples;
	int m_filterRadius;

	Timer *m_timer;

	int width, height;

	// 0: d1_min, 1: d1_max (distance from light to receiver); 2: d2_min, 3: d2_max (distance from light to occluders)
	std::vector<Vector4> dists;
	
	// 0: max_slope, 1: min_slope
	std::vector<Vector2> slopes[2];
	std::vector<bool> occluded[2];
	
	std::vector<Vector> normals;
	std::vector<Point> positions;
	std::vector<int> objId;

	std::vector<Float> filterBeta[2];

	bool m_verbose;
};

MTS_IMPLEMENT_CLASS_S(AAFIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(AAFIntegrator, "Axis-aligned filtering based integrator");
MTS_NAMESPACE_END
