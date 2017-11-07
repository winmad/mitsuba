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
#include <mitsuba/render/volume.h>
#include <mitsuba/core/statistics.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

/**
 * \brief When the following line is uncommented, the medium implementation
 * stops integrating density when it is determined that the segment has a
 * throughput of less than 'Epsilon' (see \c mitsuba/core/constants.h)
 */
#define HETVOL_EARLY_EXIT 1

/// Generate a few statistics related to the implementation?
// #define HETVOL_STATISTICS 1

#if defined(HETVOL_STATISTICS)
static StatsCounter avgNewtonIterations("Heterogeneous volume",
		"Avg. # of Newton-Bisection iterations", EAverage);
static StatsCounter avgRayMarchingStepsTransmittance("Heterogeneous volume",
		"Avg. # of ray marching steps (transmittance)", EAverage);
static StatsCounter avgRayMarchingStepsSampling("Heterogeneous volume",
		"Avg. # of ray marching steps (sampling)", EAverage);
static StatsCounter earlyExits("Heterogeneous volume",
		"Number of early exits", EPercentage);
#endif

/*!\plugin{heterogeneous}{Heterogeneous participating medium}
 * \order{2}
 * \parameters{
 *     \parameter{method}{\String}{
 *         Specifies the sampling method that is used to generate
 *         scattering events within the medium.
 *         \begin{enumerate}[(i)]
 *             \item \code{simpson}: Sampling is done by inverting a
 *             deterministic quadrature rule based on composite
 *             Simpson integration over small ray segments. Benefits
 *             from the use of good sample generators (e.g. \pluginref{ldsampler}).
 *             \item \code{woodcock}: Generate samples using
 *             Woodcock tracking. This is usually faster and
 *             always unbiased, but has the disadvantages of not benefiting
 *             from good sample generators and not providing
 *             information that is required by bidirectional
 *             rendering techniques.
 *         \end{enumerate}
 *         Default: \texttt{woodcock}
 *     }
 *     \parameter{density}{\Volume}{
 *         Volumetric data source that supplies the medium densities
 *         (in inverse scene units)
 *     }
 *     \parameter{albedo}{\Volume}{
 *         Volumetric data source that supplies the
 *         single-scattering albedo
 *     }
 *     \parameter{orientation}{\Volume}{
 *         Optional: volumetric data source that supplies the
 *         local particle orientations throughout the medium
 *     }
 *     \parameter{scale}{\Float}{
 *         Optional scale factor that will be applied to the \code{density} parameter.
 *         Provided for convenience when accomodating data based on different units,
 *         or to simply tweak the density of the medium. \default{1}
 *     }
 *     \parameter{\Unnamed}{\Phase}{
 *          A nested phase function that describes the directional
 *          scattering properties of the medium. When none is specified,
 *          the renderer will automatically use an instance of
 *          \pluginref{isotropic}.
 *     }
 * }
 *
 * \renderings{
 *     \medrendering{40}{medium_heterogeneous_density_40}
 *     \medrendering{200}{medium_heterogeneous_density_200}
 *     \medrendering{1000}{medium_heterogeneous_density_1000}
 *     \vspace{-2mm}
 *     \caption{Renderings of an index-matched medium using different scale factors (\lstref{hetvolume})}
 * }
 *
 * This plugin provides a flexible heterogeneous medium implementation, which
 * acquires its data from nested \code{volume} instances. These can be
 * constant, use a procedural function, or fetch data from disk, e.g. using a
 * memory-mapped density grid. See \secref{volumes} for details on volume data
 * sources.
 *
 * Instead of allowing separate volumes to be provided for the scattering
 * and absorption parameters \code{sigmaS} and \code{sigmaA} (as is done in
 * \pluginref{homogeneous}), this class instead takes the approach of
 * enforcing a spectrally uniform value of \code{sigmaT}, which must be
 * provided using a nested scalar-valued volume named \code{density}.
 *
 * Another nested spectrum-valued \code{albedo} volume must also be provided, which is
 * used to compute the scattering coefficient $\sigma_s$ using the expression
 * $\sigma_s = \code{scale} * \code{density} * \code{albedo}$ (i.e. 'albedo' contains the
 * single-scattering albedo of the medium.
 *
 * Optionally, one can also provide an vector-valued \code{orientation} volume,
 * which contains local particle orientation that will be passed to
 * scattering models that support this, such as a the Micro-flake or
 * Kajiya-Kay phase functions.
 *
 * \vspace{4mm}
 *
 * \begin{xml}[label=lst:hetvolume,caption=A simple heterogeneous medium backed by a grid volume]
 * <!-- Declare a heterogeneous participating medium named 'smoke' -->
 * <medium type="heterogeneous" id="smoke">
 *     <string name="method" value="simpson"/>
 *
 *     <!-- Acquire density values from an external data file -->
 *     <volume name="density" type="gridvolume">
 *         <string name="filename" value="frame_0150.vol"/>
 *     </volume>
 *
 *     <!-- The albedo is constant and set to 0.9 -->
 *     <volume name="albedo" type="constvolume">
 *         <spectrum name="value" value="0.9"/>
 *     </volume>
 *
 *     <!-- Use an isotropic phase function -->
 *     <phase type="isotropic"/>
 *
 *     <!-- Scale the density values as desired -->
 *     <float name="scale" value="200"/>
 *  </medium>
 *
 * <!-- Attach the index-matched medium to a shape in the scene -->
 * <shape type="obj">
 *     <!-- Load an OBJ file, which contains a mesh version
 *          of the axis-aligned box of the volume data file -->
 *     <string name="filename" value="bounds.obj"/>
 *
 *     <!-- Reference the medium by ID -->
 *     <ref name="interior" id="smoke"/>
 *
 *     <!-- If desired, this shape could also declare
 *          a BSDF to create an index-mismatched
 *          transition, e.g.
 *
 *     <bsdf type="dielectric"/>
 *     -->
 * </shape>
 * \end{xml}
 */
class HeterogeneousMedium : public Medium {
public:
	/// Possible integration modes
	enum EIntegrationMethod {
		/**
		 * \brief Use deterministic composite Simpson quadrature both
		 * to compute transmittances, and to sample scattering locations
		 */
		ESimpsonQuadrature = 0,

		/**
		 * \brief Use stochastic Woodcock tracking. This is potentially
		 * faster and more robust, but has the disadvantage of being
		 * incompatible with bidirectional rendering methods, which
		 * usually need to know the probability of a sample.
		 */
		EWoodcockTracking
	};

	HeterogeneousMedium(const Properties &props)
		: Medium(props) {
		m_stepSize = props.getFloat("stepSize", 0);
		m_scale = props.getFloat("scale", 1);
		if (props.hasProperty("sigmaS") || props.hasProperty("sigmaA"))
			Log(EError, "The 'sigmaS' and 'sigmaA' properties are only supported by "
				"homogeneous media. Please use nested volume instances to supply "
				"these parameters");

		if (props.hasProperty("densityMultiplier"))
			Log(EError, "The 'densityMultiplier' parameter has been deprecated and is now called 'scale'.");

		std::string method = boost::to_lower_copy(props.getString("method", "woodcock"));
		if (method == "woodcock")
			m_method = EWoodcockTracking;
		else if (method == "simpson")
			m_method = ESimpsonQuadrature;
		else
			Log(EError, "Unsupported integration method \"%s\"!", method.c_str());

		/* hetero albedo */
		m_numClusters = props.getInteger("numClusters", 1);
		m_albedoScales.resize(m_numClusters);
		for (int i = 0; i < m_numClusters; i++) {
			std::string name = formatString("albedoScale%02i", i);
			m_albedoScales[i] = props.getSpectrum(name, Spectrum(1.f));
		}
	}

	/* Unserialize from a binary data stream */
	HeterogeneousMedium(Stream *stream, InstanceManager *manager)
		: Medium(stream, manager) {
		m_method = (EIntegrationMethod) stream->readInt();
		m_scale = stream->readFloat();
		m_density = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_albedo = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_orientation = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_S1 = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_S2 = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_stepSize = stream->readFloat();

		/* hetero albedo */
		m_segmentation = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_numClusters = stream->readInt();
		m_albedoScales.resize(m_numClusters);
		for (int i = 0; i < m_numClusters; i++) {
			m_albedoScales[i] = Spectrum(stream);
		}
		configure();
	}

	/* Serialize the volume to a binary data stream */
	void serialize(Stream *stream, InstanceManager *manager) const {
		Medium::serialize(stream, manager);
		stream->writeInt(m_method);
		stream->writeFloat(m_scale);
		manager->serialize(stream, m_density.get());
		manager->serialize(stream, m_albedo.get());
		manager->serialize(stream, m_orientation.get());
		manager->serialize(stream, m_S1.get());
		manager->serialize(stream, m_S2.get());
		stream->writeFloat(m_stepSize);

		/* hetero albedo */
		manager->serialize(stream, m_segmentation.get());
		stream->writeInt(m_numClusters);
		for (int i = 0; i < m_numClusters; i++) {
			m_albedoScales[i].serialize(stream);
		}
	}

	void configure() {
		Medium::configure();
		if (m_density.get() == NULL)
			Log(EError, "No density specified!");
		if (m_albedo.get() == NULL)
			Log(EError, "No albedo specified!");

		/* hetero albedo */
		if (m_numClusters > 1 && m_segmentation.get() != NULL) {
			m_useDiffAlbedoScales = true;
		}
		else {
			m_useDiffAlbedoScales = false;
		}

		m_densityAABB = m_density->getAABB();
		m_anisotropicMedium =
			m_phaseFunction->needsDirectionallyVaryingCoefficients();

		/* Assumes that the density medium does not
		   contain values greater than one! */
		m_maxDensity = m_scale * m_density->getMaximumFloatValue();
		if (m_anisotropicMedium) {
			if (m_phaseFunction->getClass()->getName() == "SGGXPhaseFunction") {
				// hard to bound... hack by microflake
				// need to implement later!
				m_maxDensity *= m_phaseFunction->sigmaDirMax();
			}
			else {
				m_maxDensity *= m_phaseFunction->sigmaDirMax();
			}
		}
		m_invMaxDensity = 1.0f/m_maxDensity;

		if (m_stepSize == 0) {
			m_stepSize = std::min(
				m_density->getStepSize(), m_albedo->getStepSize());
			if (m_orientation != NULL)
				m_stepSize = std::min(m_stepSize,
					m_orientation->getStepSize());

			if (m_stepSize == std::numeric_limits<Float>::infinity())
				Log(EError, "Unable to infer a suitable step size for deterministic "
						"integration, please specify one manually using the 'stepSize' "
						"parameter.");
		}

		if (m_anisotropicMedium && m_orientation.get() == NULL && 
			(m_S1.get() == NULL || m_S2.get() == NULL))
			Log(EError, "Cannot use anisotropic phase function: "
				"did not specify a particle orientation field!");
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(VolumeDataSource))) {
			VolumeDataSource *volume = static_cast<VolumeDataSource *>(child);

			if (name == "albedo") {
				Assert(volume->supportsSpectrumLookups());
				m_albedo = volume;
			} else if (name == "density") {
				Assert(volume->supportsFloatLookups());
				m_density = volume;
			} else if (name == "orientation") {
				Assert(volume->supportsVectorLookups());
				m_orientation = volume;
			} else if (name == "S1") {
				Assert(volume->supportsSpectrumLookups());
				m_S1 = volume;
			} else if (name == "S2") {
				Assert(volume->supportsSpectrumLookups());
				m_S2 = volume;
			} else if (name == "segmentation") {
				Assert(volume->supportsFloatLookups());
				m_segmentation = volume;
			} else if (name == "densityTracking") {
				Assert(volume->supportsFloatLookups());
				m_densityTracking = volume;
			} else {
				Medium::addChild(name, child);
			}
		} else {
			Medium::addChild(name, child);
		}
	}

	const VolumeDataSource *getS1() const {
		return m_S1.get();
	}

	const VolumeDataSource *getS2() const {
		return m_S2.get();
	}

	/*
	 * This function uses Simpson quadrature to compute following
	 * integral:
	 *
	 *    \int_{ray.mint}^{ray.maxt} density(ray(x)) dx
	 *
	 * The integration proceeds by splitting the function into
	 * approximately \c (ray.maxt-ray.mint)/m_stepSize segments,
	 * each of which are then approximated by a quadratic polynomial.
	 * The step size must be chosen so that this approximation is
	 * valid given the behavior of the integrand.
	 *
	 * \param ray
	 *    Ray segment to be used for the integration
	 *
	 * \return
	 *    The integrated density
	 */
	Float integrateDensity(const Ray &ray) const {
		/* Determine the ray segment, along which the
		   density integration should take place */
		Float mint, maxt;
		if (!m_densityAABB.rayIntersect(ray, mint, maxt))
			return 0.0f;

		mint = std::max(mint, ray.mint);
		maxt = std::min(maxt, ray.maxt);
		Float length = maxt-mint, maxComp = 0;

		Point p = ray(mint), pLast = ray(maxt);

		/* Ignore degenerate path segments */
		for (int i=0; i<3; ++i)
			maxComp = std::max(std::max(maxComp,
				std::abs(p[i])), std::abs(pLast[i]));
		if (length < 1e-6f * maxComp)
			return 0.0f;

		/* Compute a suitable step size */
		uint32_t nSteps = (uint32_t) std::ceil(length / m_stepSize);
		nSteps += nSteps % 2;
		const Float stepSize = length/nSteps;
		const Vector increment = ray.d * stepSize;

		#if defined(HETVOL_STATISTICS)
			avgRayMarchingStepsTransmittance.incrementBase();
			earlyExits.incrementBase();
		#endif

		/* Perform lookups at the first and last node */
		Float integratedDensity = lookupDensity(p, ray.d)
			+ lookupDensity(pLast, ray.d);

		#if defined(HETVOL_EARLY_EXIT)
			const Float stopAfterDensity = -math::fastlog(Epsilon);
			const Float stopValue = stopAfterDensity*3.0f/(stepSize
					* m_scale);
		#endif

		p += increment;

		Float m = 4;
		for (uint32_t i=1; i<nSteps; ++i) {
			integratedDensity += m * lookupDensity(p, ray.d);
			m = 6 - m;

			#if defined(HETVOL_STATISTICS)
				++avgRayMarchingStepsTransmittance;
			#endif

			#if defined(HETVOL_EARLY_EXIT)
				if (integratedDensity > stopValue) {
					// Reached the threshold -- stop early
					#if defined(HETVOL_STATISTICS)
						++earlyExits;
					#endif
					return std::numeric_limits<Float>::infinity();
				}
			#endif

			Point next = p + increment;
			if (p == next) {
				Log(EWarn, "integrateDensity(): unable to make forward progress -- "
						"round-off error issues? The step size was %e, mint=%f, "
						"maxt=%f, nSteps=%i, ray=%s", stepSize, mint, maxt, nSteps,
						ray.toString().c_str());
				break;
			}
			p = next;
		}

		return integratedDensity * m_scale
			* stepSize * (1.0f / 3.0f);
	}

	/**
	 * This function uses composite Simpson quadrature to solve the
	 * following integral equation for \a t:
	 *
	 *    \int_{ray.mint}^t density(ray(x)) dx == desiredDensity
	 *
	 * The integration proceeds by splitting the function into
	 * approximately \c (ray.maxt-ray.mint)/m_stepSize segments,
	 * each of which are then approximated by a quadratic polynomial.
	 * The step size must be chosen so that this approximation is
	 * valid given the behavior of the integrand.
	 *
	 * \param ray
	 *    Ray segment to be used for the integration
	 *
	 * \param desiredDensity
	 *    Right hand side of the above equation
	 *
	 * \param integratedDensity
	 *    Contains the final integrated density. Upon success, this value
	 *    should closely match \c desiredDensity. When the equation could
	 *    \a not be solved, the parameter contains the integrated density
	 *    from \c ray.mint to \c ray.maxt (which, in this case, must be
	 *    less than \c desiredDensity).
	 *
	 * \param t
	 *    After calling this function, \c t will store the solution of the above
	 *    equation. When there is no solution, it will be set to zero.
	 *
	 * \param densityAtMinT
	 *    After calling this function, \c densityAtMinT will store the
	 *    underlying density function evaluated at \c ray(ray.mint).
	 *
	 * \param densityAtT
	 *    After calling this function, \c densityAtT will store the
	 *    underlying density function evaluated at \c ray(t). When
	 *    there is no solution, it will be set to zero.
	 *
	 * \return
	 *    When no solution can be found in [ray.mint, ray.maxt] the
	 *    function returns \c false.
	 */
	bool invertDensityIntegral(const Ray &ray, Float desiredDensity,
			Float &integratedDensity, Float &t, Float &densityAtMinT,
			Float &densityAtT) const {
		integratedDensity = densityAtMinT = densityAtT = 0.0f;

		/* Determine the ray segment, along which the
		   density integration should take place */
		Float mint, maxt;
		if (!m_densityAABB.rayIntersect(ray, mint, maxt))
			return false;
		mint = std::max(mint, ray.mint);
		maxt = std::min(maxt, ray.maxt);
		Float length = maxt - mint, maxComp = 0;
		Point p = ray(mint), pLast = ray(maxt);

		/* Ignore degenerate path segments */
		for (int i=0; i<3; ++i)
			maxComp = std::max(std::max(maxComp,
				std::abs(p[i])), std::abs(pLast[i]));
		if (length < 1e-6f * maxComp)
			return 0.0f;

		/* Compute a suitable step size (this routine samples the integrand
		   between steps, hence the factor of 2) */
		uint32_t nSteps = (uint32_t) std::ceil(length / (2*m_stepSize));
		Float stepSize = length / nSteps,
			  multiplier = (1.0f / 6.0f) * stepSize
				  * m_scale;
		Vector fullStep = ray.d * stepSize,
			   halfStep = fullStep * .5f;

		Float node1 = lookupDensity(p, ray.d);

		if (ray.mint == mint)
			densityAtMinT = node1 * m_scale;
		else
			densityAtMinT = 0.0f;

		#if defined(HETVOL_STATISTICS)
			avgRayMarchingStepsSampling.incrementBase();
		#endif

		for (uint32_t i=0; i<nSteps; ++i) {
			Float node2 = lookupDensity(p + halfStep, ray.d),
				  node3 = lookupDensity(p + fullStep, ray.d),
				  newDensity = integratedDensity + multiplier *
						(node1+node2*4+node3);
			#if defined(HETVOL_STATISTICS)
				++avgRayMarchingStepsSampling;
			#endif
			if (newDensity >= desiredDensity) {
				/* The integrated density of the last segment exceeds the desired
				   amount -- now use the Simpson quadrature expression and
				   Newton-Bisection to find the precise location of the scattering
				   event. Note that no further density queries are performed after
				   this point; instead, the density are modeled based on a
				   quadratic polynomial that is fit to the last three lookups */

				Float a = 0, b = stepSize, x = a,
					  fx = integratedDensity - desiredDensity,
					  stepSizeSqr = stepSize * stepSize,
					  temp = m_scale / stepSizeSqr;
				int it = 1;

				#if defined(HETVOL_STATISTICS)
					avgNewtonIterations.incrementBase();
				#endif
				while (true) {
					#if defined(HETVOL_STATISTICS)
						++avgNewtonIterations;
					#endif
					/* Lagrange polynomial from the Simpson quadrature */
					Float dfx = temp * (node1 * stepSizeSqr
						- (3*node1 - 4*node2 + node3)*stepSize*x
						+ 2*(node1 - 2*node2 + node3)*x*x);
					#if 0
						cout << "Iteration " << it << ":  a=" << a << ", b=" << b
							<< ", x=" << x << ", fx=" << fx << ", dfx=" << dfx << endl;
					#endif

					x -= fx/dfx;

					if (EXPECT_NOT_TAKEN(x <= a || x >= b || dfx == 0))
						x = 0.5f * (b + a);

					/* Integrated version of the above Lagrange polynomial */
					Float intval = integratedDensity + temp * (1.0f / 6.0f) * (x *
						(6*node1*stepSizeSqr - 3*(3*node1 - 4*node2 + node3)*stepSize*x
						+ 4*(node1 - 2*node2 + node3)*x*x));
					fx = intval-desiredDensity;

					if (std::abs(fx) < 1e-6f) {
						t = mint + stepSize * i + x;
						integratedDensity = intval;
						densityAtT = temp * (node1 * stepSizeSqr
							- (3*node1 - 4*node2 + node3)*stepSize*x
							+ 2*(node1 - 2*node2 + node3)*x*x);
						return true;
					} else if (++it > 30) {
						Log(EWarn, "invertDensityIntegral(): stuck in Newton-Bisection -- "
							"round-off error issues? The step size was %e, fx=%f, dfx=%f, "
							"a=%f, b=%f", stepSize, fx, dfx, a, b);
						return false;
					}

					if (fx > 0)
						b = x;
					else
						a = x;
				}
			}

			Point next = p + fullStep;
			if (p == next) {
				Log(EWarn, "invertDensityIntegral(): unable to make forward progress -- "
						"round-off error issues? The step size was %e", stepSize);
				break;
			}
			integratedDensity = newDensity;
			node1 = node3;
			p = next;
		}

		return false;
	}

	Spectrum evalTransmittance(const Ray &ray, Sampler *sampler) const {
		if (m_method == ESimpsonQuadrature || sampler == NULL) {
			return Spectrum(math::fastexp(-integrateDensity(ray)));
		} else {
			/* When Woodcock tracking is selected as the sampling method,
			   we can use this method to get a noisy (but unbiased) estimate
			   of the transmittance */
			Float mint, maxt;
			if (!m_densityAABB.rayIntersect(ray, mint, maxt))
				return Spectrum(1.0f);
			mint = std::max(mint, ray.mint);
			maxt = std::min(maxt, ray.maxt);

			#if defined(HETVOL_STATISTICS)
				avgRayMarchingStepsTransmittance.incrementBase();
			#endif
			int nSamples = 2; /// XXX make configurable
			Float result = 0;

			for (int i=0; i<nSamples; ++i) {
				Float t = mint;
				while (true) {
					t -= math::fastlog(1-sampler->next1D()) * m_invMaxDensity;
					if (t >= maxt) {
						result += 1;
						break;
					}

					Point p = ray(t);
					Float density = lookupDensity(p, ray.d) * m_scale;

					#if defined(HETVOL_STATISTICS)
						++avgRayMarchingStepsTransmittance;
					#endif

					if (density * m_invMaxDensity > sampler->next1D())
						break;
				}
			}
			return Spectrum(result/nSamples);
		}
	}

	bool sampleDistance(const Ray &ray, MediumSamplingRecord &mRec,
			Sampler *sampler) const {
		Float integratedDensity, densityAtMinT, densityAtT;
		bool success = false;
		mRec.numLobes = 1;

		if (m_method == ESimpsonQuadrature) {
			Float desiredDensity = -math::fastlog(1-sampler->next1D());
			if (invertDensityIntegral(ray, desiredDensity, integratedDensity,
					mRec.t, densityAtMinT, densityAtT)) {
				mRec.p = ray(mRec.t);
				success = true;
				Spectrum albedo;
				int clusterIndex = 0;
				if (m_useDiffAlbedoScales) {
					albedo = m_albedo->lookupSpectrum(mRec.p);
					clusterIndex = (int)m_segmentation->lookupFloat(mRec.p);
					albedo *= m_albedoScales[clusterIndex];
				}
				else {
					albedo = m_albedo->lookupSpectrum(mRec.p) * m_albedoScales[0];
				}
				mRec.sigmaS = albedo * densityAtT;
				mRec.clusterIndex = clusterIndex;
				mRec.albedoScale = m_albedoScales[clusterIndex];
				mRec.sigmaA = Spectrum(densityAtT) - mRec.sigmaS;
				mRec.orientation = m_orientation != NULL
					? m_orientation->lookupVector(mRec.p) : Vector(0.0f);

// 				if (m_orientation != NULL) {
// 					Spectrum spect = m_orientation->lookupSpectrum(mRec.p);
// 					mRec.orientation.x = spect[0];
// 					mRec.orientation.y = spect[1];
// 					mRec.orientation.z = spect[2];
// 					if (mRec.orientation.length() > 1e-6f)
// 						mRec.orientation = normalize(mRec.orientation);
// 					else
// 						mRec.orientation = Vector(0.f);
// 				}
			}

			Float expVal = math::fastexp(-integratedDensity);
			mRec.pdfFailure = expVal;
			mRec.pdfSuccess = expVal * densityAtT;
			mRec.pdfSuccessRev = expVal * densityAtMinT;
			mRec.transmittance = Spectrum(expVal);
			mRec.time = ray.time;
		} else {
			/* The following information is invalid when
			   using Woodcock-tracking */
			mRec.pdfFailure = 1.0f;
			mRec.pdfSuccess = 1.0f;
			mRec.pdfSuccessRev = 1.0f;
			mRec.transmittance = Spectrum(1.0f);
			mRec.time = ray.time;

			#if defined(HETVOL_STATISTICS)
				avgRayMarchingStepsSampling.incrementBase();
			#endif

			Float mint, maxt;
			if (!m_densityAABB.rayIntersect(ray, mint, maxt))
				return false;
			mint = std::max(mint, ray.mint);
			maxt = std::min(maxt, ray.maxt);

			Float t = mint, densityAtT = 0;
			
			Spectrum s1;
			Spectrum s2;

			while (true) {
				mRec.prevP = ray(t);
				t -= math::fastlog(1-sampler->next1D()) * m_invMaxDensity;
				if (t >= maxt)
					break;

				Point p = ray(t);

				// todo: only for experiment
				if (m_densityTracking.get() != NULL)
					densityAtT = lookupDensityTracking(p) * m_scale;
				else
					densityAtT = lookupDensity(p, ray.d, &s1, &s2) * m_scale;

				#if defined(HETVOL_STATISTICS)
					++avgRayMarchingStepsSampling;
				#endif
				if (densityAtT * m_invMaxDensity > sampler->next1D()) {
					mRec.t = t;
					mRec.p = p;
					Spectrum albedo;
					int clusterIndex = 0;
					if (m_useDiffAlbedoScales) {
						albedo = m_albedo->lookupSpectrum(mRec.p);
						clusterIndex = (int)m_segmentation->lookupFloat(mRec.p);
						albedo *= m_albedoScales[clusterIndex];
					}
					else {
						albedo = m_albedo->lookupSpectrum(mRec.p) * m_albedoScales[0];
					}
					mRec.s1[0] = s1;
					mRec.s2[0] = s2;
					mRec.pdfLobe[0] = 1;
					mRec.lobeScales[0] = 1;
					mRec.sigmaS = albedo * densityAtT;
					mRec.clusterIndex = clusterIndex;
					mRec.albedoScale = m_albedoScales[clusterIndex];
					mRec.sigmaA = Spectrum(densityAtT) - mRec.sigmaS;
					mRec.transmittance = Spectrum(densityAtT != 0.0f ? 1.0f / densityAtT : 0);
					if (!std::isfinite(mRec.transmittance[0])) // prevent rare overflow warnings
						mRec.transmittance = Spectrum(0.0f);
					mRec.orientation = m_orientation != NULL
						? m_orientation->lookupVector(p) : Vector(0.0f);

// 					if (m_orientation != NULL) {
// 						Spectrum spect = m_orientation->lookupSpectrum(mRec.p);
// 						mRec.orientation.x = spect[0];
// 						mRec.orientation.y = spect[1];
// 						mRec.orientation.z = spect[2];
// 						if (mRec.orientation.length() > 1e-6f)
// 							mRec.orientation = normalize(mRec.orientation);
// 						else
// 							mRec.orientation = Vector(0.f);
// 					}

					mRec.medium = this;
					success = true;
					break;
				}
			}
		}
		mRec.medium = this;

		return success && mRec.pdfSuccess > 0;
	}

	void eval(const Ray &ray, MediumSamplingRecord &mRec) const {
		mRec.numLobes = 1;
		Spectrum s1;
		Spectrum s2;
		if (m_method == ESimpsonQuadrature) {
			Float expVal = math::fastexp(-integrateDensity(ray));
			Float mintDensity = lookupDensity(ray(ray.mint), ray.d) * m_scale;
			Float maxtDensity = 0.0f;
			Spectrum maxtAlbedo(0.0f);
			int clusterIndex = 0;
			if (ray.maxt < std::numeric_limits<Float>::infinity()) {
				Point p = ray(ray.maxt);
				maxtDensity = lookupDensity(p, ray.d, &s1, &s2) * m_scale;
				//maxtAlbedo = m_albedo->lookupSpectrum(p);
				if (m_useDiffAlbedoScales) {
					maxtAlbedo = m_albedo->lookupSpectrum(p);
					clusterIndex = (int)m_segmentation->lookupFloat(p);
					maxtAlbedo *= m_albedoScales[clusterIndex];
				}
				else {
					maxtAlbedo = m_albedo->lookupSpectrum(p) * m_albedoScales[0];
				}
			}
			mRec.s1[0] = s1;
			mRec.s2[0] = s2;
			mRec.pdfLobe[0] = 1;
			mRec.lobeScales[0] = 1;
			mRec.transmittance = Spectrum(expVal);
			mRec.pdfFailure = expVal;
			mRec.pdfSuccess = expVal * maxtDensity;
			mRec.pdfSuccessRev = expVal * mintDensity;
			mRec.sigmaS = maxtAlbedo * maxtDensity;
			mRec.sigmaA = Spectrum(maxtDensity) - mRec.sigmaS;
			mRec.time = ray.time;
			mRec.medium = this;
			
			mRec.clusterIndex = clusterIndex;
			mRec.albedoScale = m_albedoScales[clusterIndex];
		} else {
			Log(EError, "eval(): unsupported integration method!");
		}
	}

	bool isHomogeneous() const {
		return false;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HeterogeneousMedium[" << endl
			<< "  density = " << indent(m_density.toString()) << "," << endl
			<< "  albedo = " << indent(m_albedo.toString()) << "," << endl
			<< "  orientation = " << indent(m_orientation.toString()) << "," << endl
			<< "  stepSize = " << m_stepSize << "," << endl
			<< "  scale = " << m_scale << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
	inline Float lookupDensity(const Point &p, const Vector &d, Spectrum *s1 = NULL, Spectrum *s2 = NULL) const {
		Float density = m_density->lookupFloat(p);
		if (m_anisotropicMedium && density != 0) {
			if (m_phaseFunction->getClass()->getName() == "SGGXPhaseFunction") {
				if (m_S1 == NULL) {
					Matrix3x3 D = getPhaseFunction()->getD();
					Vector orientation = m_orientation->lookupVector(p);

// 					Spectrum spect = m_orientation->lookupSpectrum(p);
// 					orientation.x = spect[0];
// 					orientation.y = spect[1];
// 					orientation.z = spect[2];
// 					if (orientation.length() > 1e-6f)
// 						orientation = normalize(orientation);
// 					else
// 						orientation = Vector(0.f);
					
					if (orientation.isZero())
						return 0;

					Vector w3 = orientation;
					Frame frame(w3);
					Matrix3x3 basis(frame.s, frame.t, w3);
					Matrix3x3 basisT;
					basis.transpose(basisT);
					Matrix3x3 S = basis * D * basisT;
					Float Sxx = S.m[0][0], Syy = S.m[1][1], Szz = S.m[2][2];
					Float Sxy = S.m[0][1], Sxz = S.m[0][2], Syz = S.m[1][2];
					Float sqrSum = Sxx * Sxx + Syy * Syy + Szz * Szz + Sxy * Sxy + Sxz * Sxz + Syz * Syz;
					//if (!(Sxx == 0 && Syy == 0 && Szz == 0 && Sxy == 0 && Sxz == 0 && Syz == 0))
					if (fabsf(sqrSum) > 1e-6f)
						density *= m_phaseFunction->sigmaDir(d, Sxx, Syy, Szz, Sxy, Sxz, Syz);
					else
						return 0;
					return density;
				}
				else {
					Spectrum S1 = m_S1->lookupSpectrum(p);
					Spectrum S2 = m_S2->lookupSpectrum(p);

					Float Sxx = S1[0], Syy = S1[1], Szz = S1[2];
					Float Sxy = S2[0], Sxz = S2[1], Syz = S2[2];

					// handle orientation transform
					Transform volumeToWorld = m_S1->getVolumeToWorld();

					Matrix3x3 Q;
					Float eig[3];

					Matrix3x3 S(Sxx, Sxy, Sxz, Sxy, Syy, Syz, Sxz, Syz, Szz);
					S.symEig(Q, eig);
					// eig[0] < eig[1] <= eig[2]
					Vector w3(Q.m[0][0], Q.m[1][0], Q.m[2][0]);
					w3 = volumeToWorld(w3);

					if (!w3.isZero()) {
						w3 = normalize(w3);

						Vector w1(Q.m[0][1], Q.m[1][1], Q.m[2][1]);
						w1 = volumeToWorld(w1);
						w1 = normalize(w1);
						Vector w2(Q.m[0][2], Q.m[1][2], Q.m[2][2]);
						w2 = volumeToWorld(w2);
						w2 = normalize(w2);

						//Frame frame(w3);
						//Matrix3x3 basis(frame.s, frame.t, w3);

						Matrix3x3 basis(w1, w2, w3);
						Matrix3x3 D(Vector(eig[1], 0, 0), Vector(0, eig[2], 0), Vector(0, 0, eig[0]));
						Matrix3x3 basisT;
						basis.transpose(basisT);
						S = basis * D * basisT;

						Sxx = S.m[0][0]; Syy = S.m[1][1]; Szz = S.m[2][2];
						Sxy = S.m[0][1]; Sxz = S.m[0][2]; Syz = S.m[1][2];

						Float sqrSum = Sxx * Sxx + Syy * Syy + Szz * Szz + Sxy * Sxy + Sxz * Sxz + Syz * Syz;
						//if (!(Sxx == 0 && Syy == 0 && Szz == 0 && Sxy == 0 && Sxz == 0 && Syz == 0))

						if (s1 && s2) {
							(*s1)[0] = Sxx; (*s1)[1] = Syy; (*s1)[2] = Szz;
							(*s2)[0] = Sxy; (*s2)[1] = Sxz; (*s2)[2] = Syz;
						}

						if (fabsf(sqrSum) > 1e-6f)
							density *= m_phaseFunction->sigmaDir(d, Sxx, Syy, Szz, Sxy, Sxz, Syz);
						else
							return 0;
						return density;
					}
					else {
						if (s1 && s2) {
							(*s1)[0] = 0; (*s1)[1] = 0; (*s1)[2] = 0;
							(*s2)[0] = 0; (*s2)[1] = 0; (*s2)[2] = 0;
						}

						return 0;
					}
				}
			}

			Vector orientation = m_orientation->lookupVector(p);
			if (!orientation.isZero())
				density *= m_phaseFunction->sigmaDir(dot(d, orientation));
			else
				return 0;
		}
		return density;
	}

	inline Float lookupDensityTracking(const Point &p) const {
		return m_densityTracking->lookupFloat(p);
	}
protected:
	EIntegrationMethod m_method;
	ref<VolumeDataSource> m_density;
	ref<VolumeDataSource> m_albedo;
	ref<VolumeDataSource> m_orientation;
	ref<VolumeDataSource> m_S1;
	ref<VolumeDataSource> m_S2;
	Float m_scale;
	bool m_anisotropicMedium;
	Float m_stepSize;
	AABB m_densityAABB;
	Float m_maxDensity;
	Float m_invMaxDensity;

	bool m_useDiffAlbedoScales;
	int m_numClusters;
	ref<VolumeDataSource> m_segmentation;
	std::vector<Spectrum> m_albedoScales;

	// only for simple experiment 
	// assume m_densityTracking share the same AABB with m_density
	// assume isotropic phase function
	ref<VolumeDataSource> m_densityTracking;
};

MTS_IMPLEMENT_CLASS_S(HeterogeneousMedium, false, Medium)
MTS_EXPORT_PLUGIN(HeterogeneousMedium, "Heterogeneous medium");
MTS_NAMESPACE_END
