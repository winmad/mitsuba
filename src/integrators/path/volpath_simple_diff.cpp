/*
	Added by Lifan Wu
	Dec 21, 2015
	
	keep tracking of derivatives with respect to albedo
*/

/*
	This file is part of Mitsuba, a physically based rendering system.

	Copyright (c) 2007-2014 by Wenzel Jakob and others.

	Mitsuba is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License Version 3
	as published by the Free Software Foundation.

	Mitsuba is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>

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
class SimpleDiffVolumetricPathTracer : public MonteCarloIntegrator {
public:
	SimpleDiffVolumetricPathTracer(const Properties &props) : MonteCarloIntegrator(props) { }

	/// Unserialize from a binary data stream
	SimpleDiffVolumetricPathTracer(Stream *stream, InstanceManager *manager)
		: MonteCarloIntegrator(stream, manager) { }

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		MediumSamplingRecord mRec;
		RayDifferential ray(r);
		Spectrum Li(0.0f);
		int index = ray.index.x + ray.index.y * width;

		bool nullChain = true, scattered = false;
		Float eta = 1.0f;

		/* Perform the first ray intersection (or ignore if the
		intersection has already been provided). */
		rRec.rayIntersect(ray);
		Spectrum throughput(1.0f);

		if (m_maxDepth == 1)
			rRec.type &= RadianceQueryRecord::EEmittedRadiance;

		/**
		* Note: the logic regarding maximum path depth may appear a bit
		* strange. This is necessary to get this integrator's output to
		* exactly match the output of other integrators under all settings
		* of this parameter.
		*/
		while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
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
				TdA[index] = throughput + TdA[index] * val;

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

						Spectrum val = value * phase->eval(
							PhaseFunctionSamplingRecord(mRec, -ray.d, dRec.d, useSGGX));
						Li += throughput * val;
						LdA[index] += TdA[index] * val;
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

				Float phaseVal = phase->sample(pRec, rRec.sampler);
				if (phaseVal == 0)
					break;
				
				throughput *= phaseVal;
				TdA[index] *= phaseVal;

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
					TdA[index] *= val;
				}

				if (!its.isValid()) {
					/* If no intersection could be found, possibly return
					attenuated radiance from a background luminaire */
					if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
						&& (!m_hideEmitters || scattered)) {
						Spectrum val = scene->evalEnvironment(ray);
						Spectrum value = throughput * val;
						if (rRec.medium) {
							Spectrum tmp = rRec.medium->evalTransmittance(ray);
							value *= tmp;
							val *= tmp;
						}

						Li += value;
						LdA[index] += TdA[index] * val;
					}
					break;
				}

				/* Possibly include emitted radiance if requested */
				if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
					&& (!m_hideEmitters || scattered)) {
					Spectrum val = its.Le(-ray.d);
					Li += throughput * val;
					LdA[index] += TdA[index] * val;
				}

				/* Include radiance from a subsurface integrator if requested */
				if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance)) {
					Spectrum val = its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);
					Li += throughput * val;
					LdA[index] += TdA[index] * val;
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
							LdA[index] += TdA[index] * val;
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
				TdA[index] *= bsdfVal;

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
				TdA[index] /= q;
			}
		}
		avgPathLength.incrementBase();
		avgPathLength += rRec.depth;
		return Li;
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int sensorResID, int samplerResID) {
		if (!MonteCarloIntegrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID))
			return false;

		Sensor *sensor = static_cast<Sensor*>(Scheduler::getInstance()->getResource(sensorResID));

		spp = sensor->getSampler()->getSampleCount();
		Log(EInfo, "spp = %d", spp);

		Vector2i filmSize = sensor->getFilm()->getSize();

		pixelNum = filmSize.x * filmSize.y;
		width = filmSize.x;
		height = filmSize.y;

		LdA = new Spectrum[pixelNum];
		TdA = new Spectrum[pixelNum];

		for (int i = 0; i < pixelNum; i++) {
			LdA[i] = Spectrum(0.f);
			TdA[i] = Spectrum(0.f);
		}

		return true;
	}

	void postprocess(const Scene *scene, RenderQueue *queue,
		const RenderJob *job, int sceneResID, int sensorResID,
		int samplerResID) {
		const Film *film = scene->getSensor()->getFilm();
		const Bitmap *bitmap = film->getImageBlock()->getBitmap();
		ref<Bitmap> img = new Bitmap(*bitmap);

		Float *data = new Float[pixelNum * 3];
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				Spectrum color(img->getPixel(Point2i(x, y)));
				for (int c = 0; c < 3; c++) {
					data[3 * (x + y * width) + c] = color[c];
				}
			}
		}
		savePfm("origin.pfm", data, width, height);

		for (int i = 0; i < pixelNum; i++) {
			LdA[i] /= (Float)spp;
			TdA[i] /= (Float)spp;
		}

		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int index = x + y * width;
				Spectrum color(TdA[index]);
				for (int c = 0; c < 3; c++) {
					data[3 * (x + y * width) + c] = color[c];
				}
			}
		}
		savePfm("TdA.pfm", data, width, height);

		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int index = x + y * width;
				Spectrum color(LdA[index]);
				for (int c = 0; c < 3; c++) {
					data[3 * (x + y * width) + c] = color[c];
				}
			}
		}
		savePfm("LdA.pfm", data, width, height);

		delete[] data;

		delete[] TdA;
		delete[] LdA;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		MonteCarloIntegrator::serialize(stream, manager);
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

	std::string toString() const {
		std::ostringstream oss;
		oss << "SimpleVolumetricPathTracer[" << endl
			<< "  maxDepth = " << m_maxDepth << "," << endl
			<< "  rrDepth = " << m_rrDepth << "," << endl
			<< "  strictNormals = " << m_strictNormals << endl
			<< "]";
		return oss.str();
	}

	int height, width, pixelNum;
	int spp;
	Spectrum *LdA;
	Spectrum *TdA;

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(SimpleDiffVolumetricPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(SimpleDiffVolumetricPathTracer, "Simple volumetric path tracer with albedo diff");
MTS_NAMESPACE_END
