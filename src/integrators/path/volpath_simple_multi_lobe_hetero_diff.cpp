/*
	Added by Lifan Wu
	May 07, 2015

	keep tracking of derivatives with respect to multiple albedo scales
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/integrator.h>

MTS_NAMESPACE_BEGIN

static StatsCounter avgPathLength("Volumetric path tracer", "Average path length", EAverage);

/*!\plugin[volpathsimple]{volpath\_simple}{Simple volumetric path tracer}
* \order{3}
* \parameters{
*     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
*         in the generated output image (where \code{-1} corresponds to $\infty$).
*	       A value of \code{1} will only render directly visible light sources.
*	       \code{2} will lead to single-bounce (direct-only) illumination,
*	       and so on. \default{\code{-1}}
*	   }
*	   \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after
*	      which the implementation will start to use the ``russian roulette''
*	      path termination criterion. \default{\code{5}}
*	   }
*     \parameter{strictNormals}{\Boolean}{Be strict about potential
*        inconsistencies involving shading normals? See
*        page~\pageref{sec:strictnormals} for details.
*        \default{no, i.e. \code{false}}
*     }
*     \parameter{hideEmitters}{\Boolean}{Hide directly visible emitters?
*        See page~\pageref{sec:hideemitters} for details.
*        \default{no, i.e. \code{false}}
*     }
* }
*
* This plugin provides a basic volumetric path tracer that can be used to
* compute approximate solutions of the radiative transfer equation. This
* particular integrator is named ``simple'' because it does not make use of
* multiple importance sampling. This results in a potentially
* faster execution time. On the other hand, it also means that this
* plugin will likely not perform well when given a scene that contains
* highly glossy materials. In this case, please use \pluginref{volpath}
* or one of the bidirectional techniques.
*
* This integrator has special support for \emph{index-matched} transmission
* events (i.e. surface scattering events that do not change the direction
* of light). As a consequence, participating media enclosed by a stencil shape (see
* \secref{shapes} for details) are rendered considerably more efficiently when this
* shape has \emph{no}\footnote{this is what signals to Mitsuba that the boundary is
* index-matched and does not interact with light in any way. Alternatively,
* the \pluginref{mask} and \pluginref{thindielectric} BSDF can be used to specify
* index-matched boundaries that involve some amount of interaction.} BSDF assigned
* to it (as compared to, say, a \pluginref{dielectric} or \pluginref{roughdielectric} BSDF).
*
* \remarks{
*    \item This integrator performs poorly when rendering
*      participating media that have a different index of refraction compared
*      to the surrounding medium.
*    \item This integrator has difficulties rendering
*      scenes that contain relatively glossy materials (\pluginref{volpath} is preferable in this case).
*    \item This integrator has poor convergence properties when rendering
*    caustics and similar effects. In this case, \pluginref{bdpt} or
*    one of the photon mappers may be preferable.
* }
*/
class SimpleMultiLobeHeteroDiffVolumetricPathTracer : public MonteCarloIntegrator {
public:
	SimpleMultiLobeHeteroDiffVolumetricPathTracer(const Properties &props) : MonteCarloIntegrator(props) {
		height = props.getInteger("height");
		width = props.getInteger("width");
		spp = props.getInteger("spp");
		numClusters = props.getInteger("numClusters", 1);
		pixelNum = height * width;
		prefix = props.getString("prefix", "");
		m_maxScattering = props.getInteger("maxScattering", -1);
		m_numLobes = props.getInteger("SGGXlobes", 1);
	}

	/// Unserialize from a binary data stream
	SimpleMultiLobeHeteroDiffVolumetricPathTracer(Stream *stream, InstanceManager *manager)
		: MonteCarloIntegrator(stream, manager) {
		height = stream->readInt();
		width = stream->readInt();
		spp = stream->readInt();
		numClusters = stream->readInt();
		prefix = stream->readString();
		m_maxScattering = stream->readInt();
		m_numLobes = stream->readInt();
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		MonteCarloIntegrator::serialize(stream, manager);
		stream->writeInt(height);
		stream->writeInt(width);
		stream->writeInt(spp);
		stream->writeInt(numClusters);
		stream->writeString(prefix);
		stream->writeInt(m_maxScattering);
		stream->writeInt(m_numLobes);
	}

	void configure() {
		pixelNum = height * width;

		LdA.resize(numClusters);
		//TdA.resize(numClusters);

		for (int i = 0; i < numClusters; i++) {
			LdA[i] = new Spectrum[pixelNum];
			//TdA[i] = new Spectrum[pixelNum];

			for (int j = 0; j < pixelNum; j++) {
				LdA[i][j] = Spectrum(0.f);
				//TdA[i][j] = Spectrum(0.f);
			}
		}

		imageSeg = new int[pixelNum];
		for (int i = 0; i < pixelNum; i++)
			imageSeg[i] = 0;

		LdW.resize(numClusters);
		for (int i = 0; i < numClusters; i++) {
			LdW[i].resize(m_numLobes);

			for (int l = 0; l < m_numLobes; l++) {
				LdW[i][l] = new Spectrum[pixelNum];

				for (int p = 0; p < pixelNum; p++) {
					LdW[i][l][p] = Spectrum(0.f);
				}
			}
		}
	}

	~SimpleMultiLobeHeteroDiffVolumetricPathTracer() {
		for (int i = 0; i < numClusters; i++) {
			if (LdA[i] != NULL)
				delete[] LdA[i];
			//if (TdA[i] != NULL)
			//	delete[] TdA[i];
		}

		if (imageSeg != NULL)
			delete[] imageSeg;

		for (int i = 0; i < numClusters; i++) {
			for (int l = 0; l < m_numLobes; l++) {
				if (LdW[i][l] != NULL)
					delete[] LdW[i][l];
			}
		}
	}

	void renderBlock(const Scene *scene,
		const Sensor *sensor, Sampler *sampler, ImageBlock *block,
		const bool &stop, const std::vector< TPoint2<uint8_t> > &points) const {

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

		std::vector<Float> cntLdA(numClusters);
		std::vector<std::vector<Float> > cntLdW;
		cntLdW.resize(numClusters);
		for (int k = 0; k < numClusters; k++) {
			cntLdW[k].resize(m_numLobes);
			for (int l = 0; l < m_numLobes; l++) {
				cntLdW[k][l] = 0.f;
			}
		}

		for (size_t i = 0; i < points.size(); ++i) {
			Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
			int index = offset.x + offset.y * width;

			if (stop)
				break;

			sampler->generate(offset);

			for (int k = 0; k < numClusters; k++) {
				cntLdA[k] = 0.f;
			}

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

				std::vector<Spectrum> oneTdA(numClusters, Spectrum(0.f));
				std::vector<Spectrum> oneLdA(numClusters, Spectrum(0.f));
				std::vector<std::vector<Spectrum> > oneTdW;
				std::vector<std::vector<Spectrum> > oneLdW;
				oneTdW.resize(numClusters);
				oneLdW.resize(numClusters);
				for (int k = 0; k < numClusters; k++) {
					oneTdW[k].resize(m_numLobes);
					oneLdW[k].resize(m_numLobes);
					for (int l = 0; l < m_numLobes; l++) {
						oneTdW[k][l] = Spectrum(0.f);
						oneLdW[k][l] = Spectrum(0.f);
					}
				}
				int albedoSegs = 0;

				spec *= Li(sensorRay, rRec, oneTdA, oneLdA, oneTdW, oneLdW, albedoSegs);

				block->put(samplePos, spec, rRec.alpha);

				for (int k = 0; k < numClusters; k++) {
					bool goodSample = true;

					for (int c = 0; c < 3; c++) {
						if (!std::isfinite(oneLdA[k][c]) || oneLdA[k][c] < 0) {
							goodSample = false;
							break;
						}
					}

					if (goodSample) {
						LdA[k][index] += oneLdA[k];
						cntLdA[k] += 1.f;
					}
				}

				for (int k = 0; k < numClusters; k++) {
					for (int l = 0; l < m_numLobes; l++) {
						bool goodSample = true;

						for (int c = 0; c < 3; c++) {
							if (!std::isfinite(oneLdW[k][l][c]) || oneLdW[k][l][c] < 0) {
								goodSample = false;
								break;
							}
						}

						if (goodSample) {
							LdW[k][l][index] += oneLdW[k][l];
							cntLdW[k][l] += 1.f;
						}
					}
				}

				imageSeg[index] |= albedoSegs;

				sampler->advance();
			}

			for (int k = 0; k < numClusters; k++) {
				if (cntLdA[k] > 0.f) {
					LdA[k][index] /= cntLdA[k];
				}
				else {
					LdA[k][index] = Spectrum(0.f);
				}

				for (int l = 0; l < m_numLobes; l++) {
					if (cntLdW[k][l] > 0.f) {
						LdW[k][l][index] /= cntLdW[k][l];
					}
					else {
						LdW[k][l][index] = Spectrum(0.f);
					}
				}
			}
		}

		Float *data = new Float[(int)points.size() * 3];

		for (int k = 0; k < numClusters; k++) {
			std::string outfile = prefix + formatString("LdA_s%02i_%03i_%03i.pfm", k, 
				block->getOffset().x, block->getOffset().y);
			for (int i = 0; i < points.size(); i++) {
				Point2i p = Point2i(points[i]);
				int localIndex = p.x + p.y * block->getWidth();
				Point2i offset = p + Vector2i(block->getOffset());
				int globalIndex = offset.x + offset.y * width;
				Spectrum color(LdA[k][globalIndex]);
				for (int c = 0; c < 3; c++) {
					data[3 * localIndex + c] = color[c];
				}
			}
			savePfm(outfile.c_str(), data, block->getWidth(), block->getHeight());
		}

		for (int k = 0; k < numClusters; k++) {
			for (int l = 0; l < m_numLobes; l++) {
				std::string outfile = prefix + formatString("LdW_s%02i_l%02i_%03i_%03i.pfm", k, l, block->getOffset().x,
					block->getOffset().y);
				for (int i = 0; i < points.size(); i++) {
					Point2i p = Point2i(points[i]);
					int localIndex = p.x + p.y * block->getWidth();
					Point2i offset = p + Vector2i(block->getOffset());
					int globalIndex = offset.x + offset.y * width;
					Spectrum color(LdW[k][l][globalIndex]);
					for (int c = 0; c < 3; c++) {
						data[3 * localIndex + c] = color[c];
					}
				}
				savePfm(outfile.c_str(), data, block->getWidth(), block->getHeight());
			}
		}

		std::string outfile = prefix + formatString("image_seg_%03i_%03i.pfm", block->getOffset().x, block->getOffset().y);
		for (int i = 0; i < points.size(); i++) {
			Point2i p = Point2i(points[i]);
			int localIndex = p.x + p.y * block->getWidth();
			Point2i offset = p + Vector2i(block->getOffset());
			int globalIndex = offset.x + offset.y * width;
			Spectrum color(imageSeg[globalIndex]);
			for (int c = 0; c < 3; c++) {
				data[3 * localIndex + c] = color[c];
			}
		}
		savePfm(outfile.c_str(), data, block->getWidth(), block->getHeight());

		/*
		outfile = formatString("TdA_%03i_%03i.pfm", block->getOffset().x, block->getOffset().y);
		for (int i = 0; i < points.size(); i++) {
		Point2i p = Point2i(points[i]);
		int localIndex = p.x + p.y * block->getWidth();
		Point2i offset = p + Vector2i(block->getOffset());
		int globalIndex = offset.x + offset.y * width;
		Spectrum color(TdA[globalIndex] / Float(spp));
		for (int c = 0; c < 3; c++) {
		data[3 * localIndex + c] = color[c];
		}
		}
		savePfm(outfile.c_str(), data, block->getWidth(), block->getHeight());
		*/
		delete[] data;
	}

	Spectrum Li(const RayDifferential &ray,
		RadianceQueryRecord &rRec) const {
		return Spectrum(0.f);
	}

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec,
		std::vector<Spectrum> &TdA, std::vector<Spectrum> &LdA, 
		std::vector<std::vector<Spectrum> > &TdW,
		std::vector<std::vector<Spectrum> > &LdW, int &albedoSegs) const {
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		MediumSamplingRecord mRec;
		RayDifferential ray(r);
		Spectrum Li(0.0f);

		bool nullChain = true, scattered = false;
		Float eta = 1.0f;

		/* Perform the first ray intersection (or ignore if the
		intersection has already been provided). */
		rRec.rayIntersect(ray);
		Spectrum throughput(1.0f);

		if (m_maxDepth == 1)
			rRec.type &= RadianceQueryRecord::EEmittedRadiance;

		int thrAlbedoSegs = 0;

		int numScattering = 0;

		/**
		* Note: the logic regarding maximum path depth may appear a bit
		* strange. This is necessary to get this integrator's output to
		* exactly match the output of other integrators under all settings
		* of this parameter.
		*/
		while ((rRec.depth <= m_maxDepth || m_maxDepth < 0) &&
			(numScattering <= m_maxScattering || m_maxScattering < 0)) {
			/* ==================================================================== */
			/*                 Radiative Transfer Equation sampling                 */
			/* ==================================================================== */
			if (rRec.medium && rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler)) {
				/* Sample the integral
				\int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
				*/
				const PhaseFunction *phase = rRec.medium->getPhaseFunction();

				Spectrum val = mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;
				throughput *= val;

				if (thrAlbedoSegs == 0)
					thrAlbedoSegs |= (1 << mRec.clusterIndex);
				for (int k = 0; k < numClusters; k++) {
					if (k == mRec.clusterIndex) {
						TdA[k] = throughput * val / mRec.albedoScale + TdA[k] * val;
					}
					else {
						TdA[k] *= val;
					}

					for (int l = 0; l < m_numLobes; l++) {
						TdW[k][l] *= val;
					}
				}

				numScattering++;

				/* ==================================================================== */
				/*                     Direct illumination sampling                     */
				/* ==================================================================== */

				/* Estimate the single scattering component if this is requested */
				if (rRec.type & RadianceQueryRecord::EDirectMediumRadiance) {
					DirectSamplingRecord dRec(mRec.p, mRec.time);
					int maxInteractions = m_maxDepth - rRec.depth - 1;

					Spectrum value = scene->sampleAttenuatedEmitterDirect(
						dRec, rRec.medium, maxInteractions,
						rRec.nextSample2D(), rRec.sampler);

					if (!value.isZero()) {
						bool useSGGX = false;
						if (phase->getClass()->getName() == "SGGXPhaseFunction")
							useSGGX = true;

						Float weightedF[MAX_SGGX_LOBES];
						Spectrum val = value * phase->eval(
							PhaseFunctionSamplingRecord(mRec, -ray.d, dRec.d, useSGGX), weightedF);
						Li += throughput * val;

						albedoSegs |= thrAlbedoSegs;
						for (int k = 0; k < numClusters; k++) {
							LdA[k] += TdA[k] * val;

							if (k == mRec.clusterIndex) {
								for (int l = 0; l < m_numLobes; l++) {
									LdW[k][l] += TdW[k][l] * val + throughput * value * weightedF[l];
								}
							}
							else {
								for (int l = 0; l < m_numLobes; l++) {
									LdW[k][l] += TdW[k][l] * val;
								}
							}
						}
					}
				}

				/* Stop if multiple scattering was not requested, or if the path gets too long */
				if ((rRec.depth + 1 >= m_maxDepth && m_maxDepth > 0) ||
					!(rRec.type & RadianceQueryRecord::EIndirectMediumRadiance))
					break;

				/* ==================================================================== */
				/*             Phase function sampling / Multiple scattering            */
				/* ==================================================================== */

				bool useSGGX = false;
				if (phase->getClass()->getName() == "SGGXPhaseFunction")
					useSGGX = true;
				PhaseFunctionSamplingRecord pRec(mRec, -ray.d, useSGGX);

				Float phasePdf;
				Float weightedF[MAX_SGGX_LOBES];
				Float phaseVal = phase->sample(pRec, phasePdf, rRec.sampler, weightedF);
				if (phaseVal == 0)
					break;

				throughput *= phaseVal;

				for (int k = 0; k < numClusters; k++) {
					TdA[k] *= phaseVal;

					if (k == mRec.clusterIndex) {
						for (int l = 0; l < m_numLobes; l++) {
							//TdW[k][l] = throughput * mRec.pdfLobe[l] + TdW[k][l] * phaseVal;
							Float dfdW = weightedF[l] / phasePdf;
							TdW[k][l] = throughput * dfdW + TdW[k][l] * phaseVal;
						}
					}
					else {
						for (int l = 0; l < m_numLobes; l++) {
							TdW[k][l] *= phaseVal;
						}
					}
				}

				/* Trace a ray in this direction */
				ray = Ray(mRec.p, pRec.wo, ray.time);
				ray.mint = 0;
				scene->rayIntersect(ray, its);
				nullChain = false;
				scattered = true;
			}
			else {
				/* Sample
				tau(x, y) * (Surface integral). This happens with probability mRec.pdfFailure
				Account for this and multiply by the proper per-color-channel transmittance.
				*/

				if (rRec.medium) {
					Spectrum val = mRec.transmittance / mRec.pdfFailure;
					throughput *= val;

					for (int k = 0; k < numClusters; k++) {
						TdA[k] *= val;
						
						for (int l = 0; l < m_numLobes; l++) {
							TdW[k][l] *= val;
						}
					}
				}

				if (!its.isValid()) {
					/* If no intersection could be found, possibly return
					attenuated radiance from a background luminaire */
					if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
						&& (!m_hideEmitters || scattered)) {
						Spectrum val = scene->evalEnvironment(ray);
						Spectrum value = throughput * val;
						if (rRec.medium) {
							Spectrum tmp = rRec.medium->evalTransmittance(ray, rRec.sampler);
							value *= tmp;
							val *= tmp;
						}

						Li += value;

						albedoSegs |= thrAlbedoSegs;
						for (int k = 0; k < numClusters; k++) {
							LdA[k] += TdA[k] * val;

							for (int l = 0; l < m_numLobes; l++) {
								LdW[k][l] += TdW[k][l] * val;
							}
						}
					}
					break;
				}

				/* Possibly include emitted radiance if requested */
				if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
					&& (!m_hideEmitters || scattered)) {
					Spectrum val = its.Le(-ray.d);
					Li += throughput * val;

					albedoSegs |= thrAlbedoSegs;
					for (int k = 0; k < numClusters; k++) {
						LdA[k] += TdA[k] * val;

						for (int l = 0; l < m_numLobes; l++) {
							LdW[k][l] += TdW[k][l] * val;
						}
					}
				}

				/* Include radiance from a subsurface integrator if requested */
				if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance)) {
					Spectrum val = its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);
					Li += throughput * val;

					albedoSegs |= thrAlbedoSegs;
					for (int k = 0; k < numClusters; k++) {
						LdA[k] += TdA[k] * val;

						for (int l = 0; l < m_numLobes; l++) {
							LdW[k][l] += TdW[k][l] * val;
						}
					}
				}

				/* Prevent light leaks due to the use of shading normals */
				Float wiDotGeoN = -dot(its.geoFrame.n, ray.d),
					wiDotShN = Frame::cosTheta(its.wi);
				if (m_strictNormals && wiDotGeoN * wiDotShN < 0)
					break;

				/* ==================================================================== */
				/*                     Direct illumination sampling                     */
				/* ==================================================================== */

				const BSDF *bsdf = its.getBSDF(ray);

				/* Estimate the direct illumination if this is requested */
				if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
					(bsdf->getType() & BSDF::ESmooth)) {
					DirectSamplingRecord dRec(its);
					int maxInteractions = m_maxDepth - rRec.depth - 1;

					Spectrum value = scene->sampleAttenuatedEmitterDirect(
						dRec, its, rRec.medium, maxInteractions,
						rRec.nextSample2D(), rRec.sampler);

					if (!value.isZero()) {
						/* Allocate a record for querying the BSDF */
						BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));
						bRec.sampler = rRec.sampler;

						Float woDotGeoN = dot(its.geoFrame.n, dRec.d);
						/* Prevent light leaks due to the use of shading normals */
						if (!m_strictNormals ||
							woDotGeoN * Frame::cosTheta(bRec.wo) > 0) {
							Spectrum val = value * bsdf->eval(bRec);
							Li += throughput * val;
							
							albedoSegs |= thrAlbedoSegs;
							for (int k = 0; k < numClusters; k++) {
								LdA[k] += TdA[k] * val;

								for (int l = 0; l < m_numLobes; l++) {
									LdW[k][l] += TdW[k][l] * val;
								}
							}
						}
					}
				}

				/* ==================================================================== */
				/*                   BSDF sampling / Multiple scattering                */
				/* ==================================================================== */

				/* Sample BSDF * cos(theta) */
				BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
				Spectrum bsdfVal = bsdf->sample(bRec, rRec.nextSample2D());
				if (bsdfVal.isZero())
					break;

				/* Recursively gather indirect illumination? */
				int recursiveType = 0;
				if ((rRec.depth + 1 < m_maxDepth || m_maxDepth < 0) &&
					(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
					recursiveType |= RadianceQueryRecord::ERadianceNoEmission;

				/* Recursively gather direct illumination? This is a bit more
				complicated by the fact that this integrator can create connection
				through index-matched medium transitions (ENull scattering events) */
				if ((rRec.depth < m_maxDepth || m_maxDepth < 0) &&
					(rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance) &&
					(bRec.sampledType & BSDF::EDelta) &&
					(!(bRec.sampledType & BSDF::ENull) || nullChain)) {
					recursiveType |= RadianceQueryRecord::EEmittedRadiance;
					nullChain = true;
				}
				else {
					nullChain &= bRec.sampledType == BSDF::ENull;
				}

				/* Potentially stop the recursion if there is nothing more to do */
				if (recursiveType == 0)
					break;
				rRec.type = recursiveType;

				/* Prevent light leaks due to the use of shading normals */
				const Vector wo = its.toWorld(bRec.wo);
				Float woDotGeoN = dot(its.geoFrame.n, wo);
				if (woDotGeoN * Frame::cosTheta(bRec.wo) <= 0 && m_strictNormals)
					break;

				/* Keep track of the throughput, medium, and relative
				refractive index along the path */
				throughput *= bsdfVal;

				for (int k = 0; k < numClusters; k++) {
					TdA[k] *= bsdfVal;

					for (int l = 0; l < m_numLobes; l++) {
						TdW[k][l] *= bsdfVal;
					}
				}

				eta *= bRec.eta;
				if (its.isMediumTransition())
					rRec.medium = its.getTargetMedium(wo);

				/* In the next iteration, trace a ray in this direction */
				ray = Ray(its.p, wo, ray.time);
				scene->rayIntersect(ray, its);
				scattered |= bRec.sampledType != BSDF::ENull;
			}

			if (rRec.depth++ >= m_rrDepth) {
				/* Russian roulette: try to keep path weights equal to one,
				while accounting for the solid angle compression at refractive
				index boundaries. Stop with at least some probability to avoid
				getting stuck (e.g. due to total internal reflection) */

				Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
				if (rRec.nextSample1D() >= q)
					break;
				throughput /= q;

				for (int k = 0; k < numClusters; k++) {
					TdA[k] /= q;

					for (int l = 0; l < m_numLobes; l++) {
						TdW[k][l] /= q;
					}
				}
			}
		}
		avgPathLength.incrementBase();
		avgPathLength += rRec.depth;

		return Li;
	}

	void savePfm(const char *fileName, float *data, int width, int height) const {
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

	std::string toString() const {
		std::ostringstream oss;
		oss << "SimpleVolumetricPathTracer[" << endl
			<< "  maxDepth = " << m_maxDepth << "," << endl
			<< "  rrDepth = " << m_rrDepth << "," << endl
			<< "  strictNormals = " << m_strictNormals << endl
			<< "]";
		return oss.str();
	}

	int numClusters;

	int height, width, pixelNum;
	int spp;
	std::vector<Spectrum*> LdA;
	//std::vector<Spectrum*> TdA;

	int *imageSeg;

	std::string prefix;

	int m_maxScattering;

	int m_numLobes;
	std::vector<std::vector<Spectrum*> > LdW;

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(SimpleMultiLobeHeteroDiffVolumetricPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(SimpleMultiLobeHeteroDiffVolumetricPathTracer, "Simple volumetric path tracer with albedo scale and lobe scale diff, for heterogeneous medium");
MTS_NAMESPACE_END
