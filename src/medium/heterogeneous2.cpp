/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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
#include <mitsuba/render/medium.h>
#include <mitsuba/render/volume2.h>
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


class HeterogeneousMediumEx : public Medium {
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

	HeterogeneousMediumEx(const Properties &props) 
		: Medium(props) {
		m_stepSize = props.getFloat("stepSize", 0);
		m_scale = props.getFloat("scale", 1);
		if (props.hasProperty("sigmaS") || props.hasProperty("sigmaA"))
			Log(EError, "The 'sigmaS' and 'sigmaA' properties are only supported by "
				"homogeneous media. Please use nested volume instances to supply "
				"these parameters");

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

		/* multi-lobe SGGX */
		m_numLobes = props.getInteger("SGGXlobes", 1);
		m_lobeScales.resize(m_numClusters);
		for (int i = 0; i < m_numClusters; i++) {
			m_lobeScales[i].resize(m_numLobes);
			for (int j = 0; j < m_numLobes; j++) {
				std::string name = formatString("lobeScale_s%02i_l%02i", i, j);
				m_lobeScales[i][j] = props.getFloat(name, 1.f);
			}
		}
		Log(EInfo, "Number of SGGX lobes = %d", m_numLobes);
	}

	/* Unserialize from a binary data stream */
	HeterogeneousMediumEx(Stream *stream, InstanceManager *manager) 
		: Medium(stream, manager) {
		m_method = (EIntegrationMethod) stream->readInt();
		m_scale = stream->readFloat();
		m_volume = static_cast<VolumeDataSourceEx *>(manager->getInstance(stream));
		m_stepSize = stream->readFloat();

		/* hetero albedo */
		m_numClusters = stream->readInt();
		m_albedoScales.resize(m_numClusters);
		for (int i = 0; i < m_numClusters; i++) {
			m_albedoScales[i] = Spectrum(stream);
		}

		/* multi-lobe SGGX */
		m_numLobes = stream->readInt();
		m_lobeScales.resize(m_numClusters);
		for (int i = 0; i < m_numClusters; i++) {
			m_lobeScales[i].resize(m_numLobes);
			for (int j = 0; j < m_numLobes; j++) {
				m_lobeScales[i][j] = stream->readFloat();
			}
		}
		configure();
	}

    virtual ~HeterogeneousMediumEx() {
    }

	/* Serialize the volume to a binary data stream */
	void serialize(Stream *stream, InstanceManager *manager) const {
		Medium::serialize(stream, manager);
		stream->writeInt(m_method);
		stream->writeFloat(m_scale);
		manager->serialize(stream, m_volume.get());
		stream->writeFloat(m_stepSize);

		/* hetero albedo */
		stream->writeInt(m_numClusters);
		for (int i = 0; i < m_numClusters; i++) {
			m_albedoScales[i].serialize(stream);
		}

		/* multi-lobe SGGX */
		stream->writeInt(m_numLobes);
		for (int i = 0; i < m_numClusters; i++) {
			for (int j = 0; j < m_numLobes; j++) {
				stream->writeFloat(m_lobeScales[i][j]);
			}
		}
	}

	void configure() {
		Medium::configure();
		if (m_volume.get() == NULL)
			Log(EError, "No volume specified!");
		m_volumeAABB = m_volume->getAABB();
        if ( !m_volumeAABB.isValid() )
            Log(EError, "Invalid volume AABB!");

		/* hetero albedo */
		if (m_numClusters > 1) {
			m_useDiffAlbedoScales = true;
		}
		else {
			m_useDiffAlbedoScales = false;
		}

		/* Assumes that the density medium does not 
		   contain values greater than one! */
		m_maxDensity = m_scale*m_volume->getMaximumFloatValueEx(0);
		if (m_volume->hasSGGXVolume())
			m_maxDensity *= m_phaseFunction->sigmaDirMax();
		m_invMaxDensity = 1.0f/m_maxDensity;

        if (m_stepSize == 0) {
            m_stepSize = m_volume->getStepSize();
			if (m_stepSize == std::numeric_limits<Float>::infinity()) 
				Log(EError, "Unable to infer a suitable step size for deterministic "
						"integration, please specify one manually using the 'stepSize' "
						"parameter.");
		}
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
        const Class *cClass = child->getClass();

		if (cClass->derivesFrom(MTS_CLASS(VolumeDataSourceEx))) {
			VolumeDataSourceEx *volume = static_cast<VolumeDataSourceEx *>(child);
            Assert(m_volume == NULL && volume->supportsBundleLookups());
            m_volume = volume;
        } else {
			Medium::addChild(name, child);
		}
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
		if (!m_volumeAABB.rayIntersect(ray, mint, maxt))
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
		if (!m_volumeAABB.rayIntersect(ray, mint, maxt))
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
			if (!m_volumeAABB.rayIntersect(ray, mint, maxt))
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
		mRec.numLobes = m_numLobes;

		if (m_method == ESimpsonQuadrature) {
			Float desiredDensity = -math::fastlog(1-sampler->next1D());
			if (invertDensityIntegral(ray, desiredDensity, integratedDensity, 
					mRec.t, densityAtMinT, densityAtT)) {
				mRec.p = ray(mRec.t);
				success = true;
				Spectrum albedo;
				Float fClusterIndex = 0.f;
				int clusterIndex = 0;

				Spectrum s1[MAX_SGGX_LOBES];
				Spectrum s2[MAX_SGGX_LOBES];
				Float pdfLobe[MAX_SGGX_LOBES];

				if (m_volume->hasOrientation())
					m_volume->lookupBundle(mRec.p, NULL, &mRec.orientation, &albedo, NULL,
						&fClusterIndex, NULL, NULL, NULL);
				else {
					m_volume->lookupBundle(mRec.p, NULL, NULL, &albedo, NULL,
						&fClusterIndex, s1, s2, pdfLobe, false);
				}

				if (m_useDiffAlbedoScales) {
					clusterIndex = (int)fClusterIndex;
					albedo *= m_albedoScales[clusterIndex];
				}
				else {
					clusterIndex = 0;
					albedo *= m_albedoScales[0];
				}

				mRec.sigmaS = albedo * densityAtT;
				mRec.clusterIndex = clusterIndex;
				mRec.albedoScale = m_albedoScales[clusterIndex];
				mRec.sigmaA = Spectrum(densityAtT) - mRec.sigmaS;

				for (int i = 0; i < m_numLobes; i++) {
					mRec.s1[i] = s1[i];
					mRec.s2[i] = s2[i];
					mRec.pdfLobe[i] = pdfLobe[i];
					mRec.lobeScales[i] = m_lobeScales[clusterIndex][i];
				}
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
			if (!m_volumeAABB.rayIntersect(ray, mint, maxt))
				return false;
			mint = std::max(mint, ray.mint);
			maxt = std::min(maxt, ray.maxt);

			Float t = mint, densityAtT = 0;

			Spectrum s1[MAX_SGGX_LOBES];
			Spectrum s2[MAX_SGGX_LOBES];
			Float pdfLobe[MAX_SGGX_LOBES];

			while (true) {
				mRec.prevP = ray(t);
				t -= math::fastlog(1-sampler->next1D()) * m_invMaxDensity;
				if (t >= maxt)
					break;

				Point p = ray(t);
				Spectrum albedo;
				Float fClusterIndex = 0.f;
				int clusterIndex = 0;

				densityAtT = lookupDensity(p, ray.d, &albedo, &fClusterIndex, s1, s2, pdfLobe) * m_scale;
			
				#if defined(HETVOL_STATISTICS)
					++avgRayMarchingStepsSampling;
				#endif
				if (densityAtT * m_invMaxDensity > sampler->next1D()) {
					mRec.t = t;
					mRec.p = p;
					
					
					if (m_volume->hasOrientation())
						m_volume->lookupBundle(p, NULL, &mRec.orientation, &albedo, NULL,
							&fClusterIndex, NULL, NULL, NULL);
// 					else
// 						m_volume->lookupBundle(p, NULL, NULL, &albedo, NULL,
// 							&fClusterIndex, &s1, &s2, &cdfLobe);
					
					if (m_useDiffAlbedoScales) {
						clusterIndex = (int)fClusterIndex;
						albedo *= m_albedoScales[clusterIndex];
					}
					else {
						clusterIndex = 0;
						albedo *= m_albedoScales[0];
					}
					
					mRec.sigmaS = albedo * densityAtT;
					mRec.clusterIndex = clusterIndex;
					mRec.albedoScale = m_albedoScales[clusterIndex];
					mRec.sigmaA = Spectrum(densityAtT) - mRec.sigmaS;

					for (int i = 0; i < m_numLobes; i++) {
						mRec.s1[i] = s1[i];
						mRec.s2[i] = s2[i];
						mRec.pdfLobe[i] = pdfLobe[i];
						mRec.lobeScales[i] = m_lobeScales[clusterIndex][i];
					}
                    // XXX - what if a single channel has a 0 intensity value?
					mRec.transmittance = mRec.sigmaS.isZero() 
						? Spectrum(0.0f) : albedo/mRec.sigmaS;
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
		mRec.numLobes = m_numLobes;

		if (m_method == ESimpsonQuadrature) {
			Spectrum s1[MAX_SGGX_LOBES];
			Spectrum s2[MAX_SGGX_LOBES];
			Float pdfLobe[MAX_SGGX_LOBES];
			
			Float expVal = math::fastexp(-integrateDensity(ray));
			Float mintDensity = lookupDensity(ray(ray.mint), ray.d) * m_scale;
			Float maxtDensity = 0.0f;
			Spectrum maxtAlbedo(0.0f);
			Float fClusterIndex = 0.f;
			int clusterIndex = 0;
			
			if (ray.maxt < std::numeric_limits<Float>::infinity()) {
				Point p = ray(ray.maxt);

				maxtDensity = lookupDensity(p, ray.d, &maxtAlbedo, &fClusterIndex, s1, s2, pdfLobe) * m_scale;

				if (m_useDiffAlbedoScales) {
					clusterIndex = (int)fClusterIndex;
					maxtAlbedo *= m_albedoScales[clusterIndex];
				}
				else {
					clusterIndex = 0;
					maxtAlbedo *= m_albedoScales[0];
				}
			}
			mRec.transmittance = Spectrum(expVal);
			mRec.pdfFailure = expVal;
			mRec.pdfSuccess = expVal * maxtDensity;
			mRec.pdfSuccessRev = expVal * mintDensity;
			mRec.sigmaS = maxtAlbedo * maxtDensity;
			mRec.sigmaA = Spectrum(maxtDensity) - mRec.sigmaS;

			for (int i = 0; i < m_numLobes; i++) {
				mRec.s1[i] = s1[i];
				mRec.s2[i] = s2[i];
				mRec.pdfLobe[i] = pdfLobe[i];
				mRec.lobeScales[i] = m_lobeScales[clusterIndex][i];
			}

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

	const VolumeDataSourceEx *getShellmap() const {
		return m_volume.get();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HeterogeneousMediumEx[" << endl
			<< "  volume = " << indent(m_volume.toString()) << "," << endl
			<< "  stepSize = " << m_stepSize << "," << endl
			<< "  densityMultiplier = " << m_scale << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()

protected:
	inline Float lookupDensity(const Point &p, const Vector &d, 
		Spectrum *albedo = NULL, Float *clusterIndex = NULL,
		Spectrum *s1 = NULL, Spectrum *s2 = NULL, Float *pdfLobe = NULL) const {
        Float density;
        Vector orientation;
		Float _clusterIndex;

		Spectrum S1[MAX_SGGX_LOBES];
		Spectrum S2[MAX_SGGX_LOBES];
		Float _pdfLobe[MAX_SGGX_LOBES];
		Float weightedPdfLobe[MAX_SGGX_LOBES];

		if (clusterIndex) *clusterIndex = 0.f;

		if (m_phaseFunction->getClass()->getName() == "SGGXPhaseFunction")  {
			if (!m_volume->hasSGGXVolume()) {
				Matrix3x3 D = getPhaseFunction()->getD();
				m_volume->lookupBundle(p, &density, &orientation, albedo, NULL,
					clusterIndex, NULL, NULL, NULL);

				if (density == 0 || orientation.isZero())
					return 0.f;

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
					return 0.f;
				return density;
			}
			else {
// 				bool lazy = false;
// 				if (m_volume->getClass()->getName() == "ShellMappedDataSourceEx")
// 					lazy = false;

				bool lazy = true;

				int clusterIdx;

				if (m_useDiffAlbedoScales) {
					m_volume->lookupBundle(p, &density, NULL, albedo, NULL,
						&_clusterIndex, S1, S2, _pdfLobe, lazy);
					clusterIdx = (int)_clusterIndex;
				}
				else {
					m_volume->lookupBundle(p, &density, NULL, albedo, NULL,
						clusterIndex, S1, S2, _pdfLobe, lazy);
					clusterIdx = 0;
				}

				if (density == 0)
					return 0;

				if (clusterIndex && m_useDiffAlbedoScales) 
					*clusterIndex = _clusterIndex;

				for (int i = 0; i < m_numLobes; i++) {
					weightedPdfLobe[i] = _pdfLobe[i] * m_lobeScales[clusterIdx][i];
				}

				if (!lazy) {
					density *= m_phaseFunction->sigmaDir(d, S1, S2, weightedPdfLobe, m_numLobes);
				}
				else {
					Vector w1[MAX_SGGX_LOBES];
					Vector w2[MAX_SGGX_LOBES];
					Vector w3[MAX_SGGX_LOBES];
					Vector sigmaSqr[MAX_SGGX_LOBES];

					m_volume->lookupSGGXFrame(p, w1, w2, w3, sigmaSqr);

					Matrix3x3 basisT;
					Matrix3x3 S;

					for (int i = 0; i < m_numLobes; i++) {
						if (_pdfLobe[i] > 0 && !w3[i].isZero()) {
							Matrix3x3 basis(w1[i], w2[i], w3[i]);
							Matrix3x3 D(Vector(sigmaSqr[i][0], 0, 0), 
								Vector(0, sigmaSqr[i][1], 0), 
								Vector(0, 0, sigmaSqr[i][2]));

							basis.transpose(basisT);
							S = basis * D * basisT;

							S1[i][0] = S.m[0][0]; S1[i][1] = S.m[1][1]; S1[i][2] = S.m[2][2];
							S2[i][0] = S.m[0][1]; S2[i][1] = S.m[0][2]; S2[i][2] = S.m[1][2];
						}
						else {
							S1[i] = Spectrum(0.f);
							S2[i] = Spectrum(0.f);
						}
						
					}

					density *= m_phaseFunction->sigmaDir(d, S1, S2, weightedPdfLobe, m_numLobes);
				}

				if (s1 && s2 && pdfLobe) {
					for (int i = 0; i < m_numLobes; i++) {
						s1[i] = S1[i];
						s2[i] = S2[i];
						pdfLobe[i] = _pdfLobe[i];
					}
				}

				return density;
			}
		}
		else {
			m_volume->lookupBundle(p, &density, &orientation, albedo, NULL,
				clusterIndex, NULL, NULL, NULL);
			if (density != 0 && !orientation.isZero())
				return density*m_phaseFunction->sigmaDir(dot(d, orientation));
			else
				return 0.0f;
		}
	}

protected:
	EIntegrationMethod m_method;

    ref<VolumeDataSourceEx> m_volume;
	Float m_scale;
	Float m_stepSize;
	AABB m_volumeAABB;
	Float m_maxDensity;
	Float m_invMaxDensity;

	bool m_useDiffAlbedoScales;
	int m_numClusters;
	std::vector<Spectrum> m_albedoScales;

	int m_numLobes;
	std::vector<std::vector<Float> > m_lobeScales;
};

MTS_IMPLEMENT_CLASS_S(HeterogeneousMediumEx, false, Medium)
MTS_EXPORT_PLUGIN(HeterogeneousMediumEx, "Heterogeneous medium 2");
MTS_NAMESPACE_END
