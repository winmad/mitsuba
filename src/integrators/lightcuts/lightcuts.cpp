/*
	Added by Lifan Wu
	Oct 18, 2015
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/render/lightcutter.h>
#include <mitsuba/render/vpl.h>

MTS_NAMESPACE_BEGIN

#define DIRECT 0

class LightcutsIntegrator : public SamplingIntegrator {
public:
	LightcutsIntegrator(const Properties &props) : SamplingIntegrator(props) {
		/* Number of virtual lights */
		m_vplSamples = props.getSize("vplSamples", 1000);
		m_maxCutSize = props.getSize("maxCutSize", 300);
		m_maxErrorRatio = props.getFloat("maxErrorRatio", 0.01f);
		m_gLimit = props.getFloat("gLimit", 1e-3f);
		m_dLimit = props.getFloat("dLimit", 0.1f);
		m_useCosBound = props.getBoolean("useCosBound", false);
		/* Number of shading samples -- this parameter is a shorthand notation
		to set both 'emitterSamples' and 'bsdfSamples' at the same time*/
		size_t shadingSamples = props.getSize("shadingSamples", 1);

		/* Number of samples to take using the emitter sampling technique */
		m_emitterSamples = props.getSize("emitterSamples", shadingSamples);
		/* Number of samples to take using the BSDF sampling technique */
		m_bsdfSamples = props.getSize("bsdfSamples", shadingSamples);
		/* Be strict about potential inconsistencies involving shading normals? */
		m_strictNormals = props.getBoolean("strictNormals", false);
		/* When this flag is set to true, contributions from directly
		* visible emitters will not be included in the rendered image */
		m_hideEmitters = props.getBoolean("hideEmitters", false);
		Assert(m_emitterSamples + m_bsdfSamples > 0);
	}

	/// Unserialize from a binary data stream
	LightcutsIntegrator(Stream *stream, InstanceManager *manager)
		: SamplingIntegrator(stream, manager) {
		m_emitterSamples = stream->readSize();
		m_bsdfSamples = stream->readSize();
		m_strictNormals = stream->readBool();
		m_hideEmitters = stream->readBool();
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		SamplingIntegrator::serialize(stream, manager);
		stream->writeSize(m_emitterSamples);
		stream->writeSize(m_bsdfSamples);
		stream->writeBool(m_strictNormals);
		stream->writeBool(m_hideEmitters);
	}

	void configure() {
		SamplingIntegrator::configure();

		size_t sum = m_emitterSamples + m_bsdfSamples;
		m_weightBSDF = 1 / (Float)m_bsdfSamples;
		m_weightLum = 1 / (Float)m_emitterSamples;
		m_fracBSDF = m_bsdfSamples / (Float)sum;
		m_fracLum = m_emitterSamples / (Float)sum;
	}

	void configureSampler(const Scene *scene, Sampler *sampler) {
		SamplingIntegrator::configureSampler(scene, sampler);
		if (m_emitterSamples > 1)
			sampler->request2DArray(m_emitterSamples);
		if (m_bsdfSamples > 1)
			sampler->request2DArray(m_bsdfSamples);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue,
		const RenderJob *job, int sceneResID, int sensorResID,
		int samplerResID) {
		Integrator::preprocess(scene, queue, job, sceneResID,
			sensorResID, samplerResID);

		if (!(scene->getSensor()->getType() & Sensor::EProjectiveCamera))
			Log(EError, "The VPL integrator requires a projective camera "
			"(e.g. perspective/thinlens/orthographic/telecentric)!");

		m_random = new Random();
		m_ratio = scene->getBSphere().radius / 16.f;

		m_pointLightTree.init(m_random, 0.f);
		m_directionalLightTree.init(m_random, 0.f);
		m_surfaceLightTree.init(m_random, m_ratio);

		std::deque<VPL> _vpls;
		m_surfaceVPLs.clear();
		Float normalization = 1.f / generateVPLs(scene, m_random, 0, m_vplSamples, 5, true, _vpls);
		for (int i = 0; i < _vpls.size(); i++) {
			_vpls[i].P *= normalization;
			if (_vpls[i].type == ESurfaceVPL) {
				m_surfaceVPLs.push_back(_vpls[i]);
			}
		}
		
		Log(EInfo, "Building light tree begin");
		m_surfaceLightTree.build(m_surfaceVPLs);
		Log(EInfo, "Finish building light tree");

		//m_lightcutter.init(&m_pointLightTree, &m_directionalLightTree, &m_surfaceLightTree,
		//	1e-4f, 1e-4f);

		m_lightcutter = new Lightcutter(&m_pointLightTree, &m_directionalLightTree, &m_surfaceLightTree,
			m_gLimit, m_dLimit, 1e-4f, m_useCosBound);

		return true;
	}
        
        bool render(Scene *scene,
		RenderQueue *queue, const RenderJob *job,
		int sceneResID, int sensorResID, int samplerResID)
        {
            SamplingIntegrator::render(scene, queue, job, sceneResID, sensorResID, samplerResID);
        }

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		RayDifferential ray(r);
		Spectrum Li(0.0f);
		Point2 sample;

		/* Perform the first ray intersection (or ignore if the
		intersection has already been provided). */
		if (!rRec.rayIntersect(ray)) {
			/* If no intersection could be found, possibly return
			radiance from a background emitter */
			if (rRec.type & RadianceQueryRecord::EEmittedRadiance && !m_hideEmitters)
				return scene->evalEnvironment(ray);
			else
				return Spectrum(0.0f);
		}

		/* Possibly include emitted radiance if requested */
		if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance) && !m_hideEmitters)
			Li += its.Le(-ray.d);

		/* Include radiance from a subsurface scattering model if requested */
		if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
			Li += its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

		const BSDF *bsdf = its.getBSDF(ray);

		if (!(rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)
			|| (m_strictNormals && dot(ray.d, its.geoFrame.n)
			* Frame::cosTheta(its.wi) >= 0)) {
			/* Only render the direct illumination component if
			*
			* 1. It was requested
			* 2. The surface has an associated BSDF (i.e. it isn't an index-
			*    matched medium transition -- this is not supported by 'direct')
			* 3. If 'strictNormals'=true, when the geometric and shading
			*    normals classify the incident direction to the same side
			*/
			return Li;
		}
		
		/* ==================================================================== */
		/*                          Emitter sampling                          */
		/* ==================================================================== */
		bool adaptiveQuery = (rRec.extra & RadianceQueryRecord::EAdaptiveQuery);

		/* Figure out how many BSDF and direct illumination samples to
		generate, and where the random numbers should come from */
		Point2 *sampleArray;
		size_t numDirectSamples = m_emitterSamples,
			numBSDFSamples = m_bsdfSamples;
		Float fracLum = m_fracLum, fracBSDF = m_fracBSDF,
			weightLum = m_weightLum, weightBSDF = m_weightBSDF;

		if (rRec.depth > 1 || adaptiveQuery) {
			/* This integrator is used recursively by another integrator.
			Be less accurate as this sample will not directly be observed. */
			numBSDFSamples = numDirectSamples = 1;
			fracLum = fracBSDF = .5f;
			weightLum = weightBSDF = 1.0f;
		}

		if (numDirectSamples > 1) {
			sampleArray = rRec.sampler->next2DArray(numDirectSamples);
		}
		else {
			sample = rRec.nextSample2D(); sampleArray = &sample;
		}

		DirectSamplingRecord dRec(its);
		if (bsdf->getType() & BSDF::ESmooth) {
			/* Only use direct illumination sampling when the surface's
			BSDF has smooth (i.e. non-Dirac delta) component */
			for (size_t i = 0; i<numDirectSamples; ++i) {
				/* Estimate the direct illumination if this is requested */
				Spectrum value = scene->sampleEmitterDirect(dRec, sampleArray[i]);
                                
                                Intersection occ_its;
                                Ray ray(dRec.ref, dRec.d, Epsilon,
					dRec.dist*(1-ShadowEpsilon), dRec.time);
                                bool hit = scene->rayIntersect(ray, occ_its);
                                if(hit && occ_its.t < its.t)
                                    return Spectrum((its.t/(its.t-occ_its.t)-1));
                                else
                                    return Spectrum(0.0);
                                
				if (!value.isZero()) {
					const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

					/* Allocate a record for querying the BSDF */
					BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));

					/* Evaluate BSDF * cos(theta) */
					const Spectrum bsdfVal = bsdf->eval(bRec);

					if (!bsdfVal.isZero() && (!m_strictNormals
						|| dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {
						/* Calculate prob. of sampling that direction using BSDF sampling */
						Float bsdfPdf = emitter->isOnSurface() ? bsdf->pdf(bRec) : 0;

						/* Weight using the power heuristic */
						const Float weight = miWeight(dRec.pdf * fracLum,
							bsdfPdf * fracBSDF) * weightLum;

#if DIRECT
						Li += value * bsdfVal * weight;
#endif
					}
				}
			}
		}

		/* ==================================================================== */
		/*                            BSDF sampling                             */
		/* ==================================================================== */

		if (numBSDFSamples > 1) {
			sampleArray = rRec.sampler->next2DArray(numBSDFSamples);
		}
		else {
			sample = rRec.nextSample2D(); sampleArray = &sample;
		}

		Intersection bsdfIts;
		for (size_t i = 0; i<numBSDFSamples; ++i) {
			/* Sample BSDF * cos(theta) and also request the local density */
			Float bsdfPdf;

			BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
			Spectrum bsdfVal = bsdf->sample(bRec, bsdfPdf, sampleArray[i]);
			if (bsdfVal.isZero())
				continue;

			/* Prevent light leaks due to the use of shading normals */
			const Vector wo = its.toWorld(bRec.wo);
			Float woDotGeoN = dot(its.geoFrame.n, wo);
			if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
				continue;

			/* Trace a ray in this direction */
			Ray bsdfRay(its.p, wo, ray.time);

			Spectrum value;
			if (scene->rayIntersect(bsdfRay, bsdfIts)) {
				/* Intersected something - check if it was an emitter */
				if (!bsdfIts.isEmitter())
					continue;

				value = bsdfIts.Le(-bsdfRay.d);
				dRec.setQuery(bsdfRay, bsdfIts);
			}
			else {
				/* Intersected nothing -- perhaps there is an environment map? */
				const Emitter *env = scene->getEnvironmentEmitter();

				if (!env || (m_hideEmitters && bRec.sampledType == BSDF::ENull))
					continue;

				value = env->evalEnvironment(RayDifferential(bsdfRay));
				if (!env->fillDirectSamplingRecord(dRec, bsdfRay))
					continue;
			}

			/* Compute the prob. of generating that direction using the
			implemented direct illumination sampling technique */
			const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
				scene->pdfEmitterDirect(dRec) : 0;

			/* Weight using the power heuristic */
			const Float weight = miWeight(bsdfPdf * fracBSDF,
				lumPdf * fracLum) * weightBSDF;

#if DIRECT
			Li += value * bsdfVal * weight;
#endif
		}

		/* Lightcuts */
		Spectrum indirResult(0.f);

		indirResult += m_lightcutter->evalLightcut(ray, rRec, m_random, m_maxCutSize, m_maxErrorRatio);
		
		/*
		Vector wi = -ray.d;
		for (int i = 0; i < m_surfaceLightTree.nextFreeNode; i++) {
			if (m_surfaceLightTree.m_nodes[i].isLeaf()) {
				const SurfaceLightNode *node = &m_surfaceLightTree.m_nodes[i];
				indirResult += m_lightcutter.evalNodeIllumination(node, rRec.scene, m_random, wi, its, bsdf);
			}
		}
		*/
		//Log(EInfo, "(%.6f, %.6f, %.6f)", indirResult[0], indirResult[1], indirResult[2]);
		Li += indirResult;

		return Li;
	}
	
	void postprocess(const Scene *scene, RenderQueue *queue,
		const RenderJob *job, int sceneResID, int sensorResID,
		int samplerResID) {
		m_lightcutter->avgCutSize /= m_lightcutter->count;
		Log(EInfo, "avgCutSize = %.6f", m_lightcutter->avgCutSize);
	}

	inline Float miWeight(Float pdfA, Float pdfB) const {
		pdfA *= pdfA; pdfB *= pdfB;
		return pdfA / (pdfA + pdfB);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "LightcutsIntegrator[" << endl
			<< "  emitterSamples = " << m_emitterSamples << "," << endl
			<< "  bsdfSamples = " << m_bsdfSamples << "," << endl
			<< "  strictNormals = " << m_strictNormals << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
public:
	size_t m_emitterSamples;
	size_t m_bsdfSamples;
	size_t m_vplSamples;
	Float m_gLimit;
	Float m_dLimit;
	int m_maxCutSize;
	Float m_maxErrorRatio;
	Float m_fracBSDF, m_fracLum;
	Float m_weightBSDF, m_weightLum;
	bool m_strictNormals;
	bool m_hideEmitters;
	bool m_useCosBound;

	std::vector<VPL> m_surfaceVPLs;

	LightTree<PointLightNode> m_pointLightTree;
	LightTree<DirectionalLightNode> m_directionalLightTree;
	LightTree<SurfaceLightNode> m_surfaceLightTree;

	Lightcutter *m_lightcutter;

	Random *m_random;
	Float m_ratio;

	void testPointLightTree() {
		std::vector<VPL> vpls;
		VPL l;
		l.its.p = Point(0, 0, 0);
		l.P = Spectrum(100);
		vpls.push_back(l);
		
		l.its.p = Point(1, 0, 0);
		l.P = Spectrum(1);
		vpls.push_back(l);
		
		l.its.p = Point(2, 0, 0);
		l.P = Spectrum(1);
		vpls.push_back(l);

		l.its.p = Point(-1, 0, 0);
		l.P = Spectrum(1);
		vpls.push_back(l);

		l.its.p = Point(-2, 0, 0);
		l.P = Spectrum(1);
		vpls.push_back(l);

		m_pointLightTree.build(vpls);
		
		std::queue<PointLightNode*> q;
		q.push(m_pointLightTree.root);
		while (!q.empty()) {
			PointLightNode *node = q.front();
			q.pop();
			VPL *light = node->light;
			if (!node->isLeaf()) {
				q.push(node->left);
				q.push(node->right);
			}
			Log(EInfo, "======= light cluster =======\n");
			Log(EInfo, "pos = (%.6f, %.6f, %.6f)\n", light->its.p.x, light->its.p.y, light->its.p.z);
			Log(EInfo, "intensity = (%.6f, %.6f, %.6f)\n", node->P[0], node->P[1], node->P[2]);
		}
	}
};

MTS_IMPLEMENT_CLASS_S(LightcutsIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(LightcutsIntegrator, "Lightcuts integrator");
MTS_NAMESPACE_END