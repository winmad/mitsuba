#include <mitsuba/render/scene.h>
#include <mitsuba/render/volume2.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/fresolver.h>
#include <boost/algorithm/string.hpp>
#include "../volume/tetra2.h"


MTS_NAMESPACE_BEGIN


#define HETVOL_EARLY_EXIT 1

/// Generate a few statistics related to the implementation?
//#define HETVOL_STATISTICS 1

#ifdef HETVOL_STATISTICS
static StatsCounter avgNewtonIterations("Heterogeneous volume", 
        "Avg. # of Newton-Bisection iterations", EAverage);
static StatsCounter avgRayMarchingStepsTransmittance("Heterogeneous volume", 
        "Avg. # of ray marching steps (transmittance)", EAverage);
static StatsCounter avgRayMarchingStepsSampling("Heterogeneous volume", 
        "Avg. # of ray marching steps (sampling)", EAverage);
static StatsCounter earlyExits("Heterogeneous volume", 
        "Number of early exits", EPercentage);
#endif


class HeterogeneousMedium3 : public Medium
{
public:
    HeterogeneousMedium3(const Properties &props) : Medium(props), m_ready(false)
    {
        m_scale = props.getFloat("scale", 1);
        if (props.hasProperty("sigmaS") || props.hasProperty("sigmaA"))
            Log(EError, "The 'sigmaS' and 'sigmaA' properties are only supported by "
                "homogeneous media. Please use nested volume instances to supply "
                "these parameters");

        if (props.hasProperty("densityMultiplier"))
            Log(EError, "The 'densityMultiplier' parameter has been deprecated and is now called 'scale'.");

        m_volumeToWorld = props.getTransform("toWorld", Transform());
        m_shellfile = props.getString("shellfile");
        fs::path resolved = Thread::getThread()->getFileResolver()->resolve(m_shellfile);
        if ( !m_shell.load(resolved.string().c_str()) )
            Log(EError, "Failed to load the shell file!");
        Log(EInfo, "Building a BVH for the shell mesh ..");
        m_shell.configure(m_volumeToWorld);
        Log(EInfo, "Finished with %u tetrahedra (tree depth: %u).",
            m_shell.getTetrahedronCount(), m_shell.getTreeDepth());

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

    
    HeterogeneousMedium3(Stream *stream, InstanceManager *manager) : Medium(stream, manager), m_ready(false)
    {
        m_scale = stream->readFloat();
        m_volume = static_cast<VolumeDataSourceEx *>(manager->getInstance(stream));

        m_volumeToWorld = Transform(stream);
        m_shellfile = stream->readString();
        if ( !m_shell.load(m_shellfile.c_str()) )
            Log(EError, "Failed to load the shell file!");
        Log(EInfo, "Building a BVH for the shell mesh ..");
        m_shell.configure(m_volumeToWorld);
        Log(EInfo, "Finished with %u tetrahedra (tree depth: %u).",
            m_shell.getTetrahedronCount(), m_shell.getTreeDepth());

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


    virtual ~HeterogeneousMedium3() {}


    void serialize(Stream *stream, InstanceManager *manager) const
    {
        Medium::serialize(stream, manager);
        stream->writeFloat(m_scale);
        manager->serialize(stream, m_volume.get());

        m_volumeToWorld.serialize(stream);
        stream->writeString(m_shellfile);

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


    void configure()
    {
        if ( !m_ready ) {
            Medium::configure();
            if ( m_volume.get() == NULL )
                Log(EError, "No volume specified!");
            m_volumeAABB = m_volume->getAABB();
            if ( !m_volumeAABB.isValid() )
                Log(EError, "Invalid volume AABB!");
            m_anisotropicMedium =
                m_phaseFunction->needsDirectionallyVaryingCoefficients();

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
            if ( m_anisotropicMedium )
                m_maxDensity *= m_phaseFunction->sigmaDirMax();
            m_invMaxDensity = 1.0f/m_maxDensity;
            Log(EInfo, "Medium Mode: %s", m_anisotropicMedium ? "Anisotropic" : "Isotropic");

            m_worldToVolume = m_volumeToWorld.inverse();
            m_textureToData = Transform::translate(Vector(m_volumeAABB.min))*Transform::scale(m_volumeAABB.getExtents());

            if ( m_anisotropicMedium && !m_volume->supportsVectorLookups() && !m_volume->hasSGGXVolume() )
                Log(EError, "Cannot use anisotropic phase function: "
                "did not specify a particle orientation field!");

            m_ready = true;
        }
    }


    void addChild(const std::string &name, ConfigurableObject *child)
    {
        const Class *cClass = child->getClass();

        if (cClass->derivesFrom(MTS_CLASS(VolumeDataSourceEx))) {
            VolumeDataSourceEx *volume = static_cast<VolumeDataSourceEx *>(child);
            Assert(m_volume == NULL && volume->supportsBundleLookups());
            m_volume = volume;
        } else {
            Medium::addChild(name, child);
        }
    }


    Float integrateDensity(const Ray &ray) const
    {
        Log(EError, "HeterogeneousMedium::integrateDensity() has not been implemented!");
        return 0.0f;
    }


    bool invertDensityIntegral(const Ray &ray, Float desiredDensity,
            Float &integratedDensity, Float &t, Float &densityAtMinT,
            Float &densityAtT) const
    {
        Log(EError, "HeterogeneousMedium::invertDensityIntegral() has not been implemented!");
        return false;
    }


    Spectrum evalTransmittance(const Ray &ray, Sampler *sampler) const
    {
        if (sampler == NULL) {
            Log(EError, "HeterogeneousMedium::evalTransmittance() needs a sampler!");
            return Spectrum(0.0f);
        } else {
            /* When Woodcock tracking is selected as the sampling method,
               we can use this method to get a noisy (but unbiased) estimate
               of the transmittance */

            Float mint, maxt;
            Ray _ray;

            uint32_t id0, f0, f1;
            Float t0, t1;
            if ( m_shell.lookupPoint(ray(ray.mint), id0) ) {
                mint = ray.mint; maxt = ray.maxt;
            }
            else {
                Intersection its;
                _ray = Ray(ray, ray.mint, std::numeric_limits<Float>::infinity());
                if ( !m_shell.m_btree->rayIntersect(_ray, its) || its.t > ray.maxt )
                    return Spectrum(1.0f);
                mint = its.t;
                if ( !m_shell.lookupPoint(ray(mint + Epsilon), id0) ) return Spectrum(1.0f);
                maxt = ray.maxt;
            }

            _ray = Ray(ray, -std::numeric_limits<Float>::infinity(), std::numeric_limits<Float>::infinity());
            if ( !m_shell.m_tetra[id0].rayIntersect(m_shell.m_vtxPosition, _ray, t0, f0, t1, f1) ) {
                //Log(EWarn, "badness");
                return Spectrum(1.0f);
            }
            
            #if defined(HETVOL_STATISTICS)
                avgRayMarchingStepsTransmittance.incrementBase();
            #endif
            int nSamples = 2; /// XXX make configurable
            Float result = 0;

			Float fClusterIndex = 0.f;
			int clusterIndex = 0;

			Spectrum s1[MAX_SGGX_LOBES];
			Spectrum s2[MAX_SGGX_LOBES];
			Float pdfLobe[MAX_SGGX_LOBES];

            for (int i=0; i<nSamples; ++i) {
                Float t = mint;
                int id = static_cast<int>(id0);
                Float minT = t0, maxT = t1;
                uint32_t minF = f0, maxF = f1;

                while (true) {
                    t -= math::fastlog(1-sampler->next1D()) * m_invMaxDensity;
                    if ( t > maxt ) {
                        result += 1.0f; break;
                    }

                    while ( t > maxT ) {
                        id = m_shell.m_link[4*id + maxF];
                        if ( id < 0 ) break;
                        if ( !m_shell.m_tetra[id].rayIntersect(m_shell.m_vtxPosition, _ray, minT, minF, maxT, maxF) ) {
                            id = -1; break;
                        }
                    }

                    if ( id < 0 ) {
                        result += 1.0f; break;
                    }
            
                    Point p = ray(t);
                    Float density = 0.0f;
                    Point4 bb;
                    if ( m_shell.m_tetra[id].inside(m_shell.m_vtxPosition, p, bb) ) {
                        const uint32_t *tmp = m_shell.m_tetra[id].idx;

                        Point q = m_textureToData.transformAffine(
                            m_shell.m_vtxTexcoord[tmp[0]]*bb.x +
                            m_shell.m_vtxTexcoord[tmp[1]]*bb.y +
                            m_shell.m_vtxTexcoord[tmp[2]]*bb.z +
                            m_shell.m_vtxTexcoord[tmp[3]]*bb.w
                        );
                        Normal norm = m_shell.m_vtxNormal[tmp[0]]*bb.x +
                            m_shell.m_vtxNormal[tmp[1]]*bb.y +
                            m_shell.m_vtxNormal[tmp[2]]*bb.z +
                            m_shell.m_vtxNormal[tmp[3]]*bb.w;
                        Vector dpdu = m_shell.m_vtxTangent[tmp[0]].dpdu*bb.x +
                            m_shell.m_vtxTangent[tmp[1]].dpdu*bb.y +
                            m_shell.m_vtxTangent[tmp[2]].dpdu*bb.z +
                            m_shell.m_vtxTangent[tmp[3]].dpdu*bb.w;
                        Vector dpdv = m_shell.m_vtxTangent[tmp[0]].dpdv*bb.x +
                            m_shell.m_vtxTangent[tmp[1]].dpdv*bb.y +
                            m_shell.m_vtxTangent[tmp[2]].dpdv*bb.z +
                            m_shell.m_vtxTangent[tmp[3]].dpdv*bb.w;

                        Vector orientation;
						
                        if ( m_anisotropicMedium ) {
							Frame tangFrame;
							tangFrame.s = dpdu; tangFrame.t = dpdv; tangFrame.n = norm;

							density = lookupDensity(q, ray.d, tangFrame, NULL, NULL, &fClusterIndex,
								s1, s2, pdfLobe) * m_scale;
		
							/*
                            m_volume->lookupBundle(q, &density, &orientation, NULL, NULL,
								NULL, NULL, NULL);
                            orientation = orientation.x*dpdu + orientation.y*dpdv + orientation.z*norm;
                            if ( density > 0.0f && !orientation.isZero() ) {
                                orientation = normalize(orientation);
                                density *= m_phaseFunction->sigmaDir(dot(ray.d, orientation))*m_scale;
                            }
							*/
                        }
                        else
                            density = m_volume->lookupFloat(q) * m_scale;
                    }
            
                    #ifdef HETVOL_STATISTICS
                        ++avgRayMarchingStepsTransmittance;
                    #endif

                    if (density * m_invMaxDensity > sampler->next1D()) 
                        break;
                }
            }
            return Spectrum(result/nSamples);
        }
    }


    bool sampleDistance(const Ray &ray, MediumSamplingRecord &mRec,	Sampler *sampler) const
    {
        bool success = false;
		mRec.numLobes = m_numLobes;

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
        Ray _ray;

        uint32_t id0, minF, maxF;
        Float minT, maxT;
        if ( m_shell.lookupPoint(ray(ray.mint), id0) ) {
            mint = ray.mint; maxt = ray.maxt;
        }
        else {
            Intersection its;
            _ray = Ray(ray, ray.mint, std::numeric_limits<Float>::infinity());
            if ( !m_shell.m_btree->rayIntersect(_ray, its) || its.t > ray.maxt )
                return false;
            mint = its.t;
            if ( !m_shell.lookupPoint(ray(mint + Epsilon), id0) ) return false;
            maxt = ray.maxt;
        }

        _ray = Ray(ray, -std::numeric_limits<Float>::infinity(), std::numeric_limits<Float>::infinity());
        if ( !m_shell.m_tetra[id0].rayIntersect(m_shell.m_vtxPosition, _ray, minT, minF, maxT, maxF) ) {
            //Log(EWarn, "badness");
            return false;
        }

        Float t = mint;
        int id = static_cast<int>(id0);

		Spectrum s1[MAX_SGGX_LOBES];
		Spectrum s2[MAX_SGGX_LOBES];
		Float pdfLobe[MAX_SGGX_LOBES];

        while (true) {
			mRec.prevP = ray(t);
            t -= math::fastlog(1 - sampler->next1D())*m_invMaxDensity;
            if ( t >= maxt )
                break;

            while ( t > maxT ) {
                id = m_shell.m_link[4*id + maxF];
                if ( id < 0 ) break;
                if ( !m_shell.m_tetra[id].rayIntersect(m_shell.m_vtxPosition, _ray, minT, minF, maxT, maxF) ) {
                    id = -1; break;
                }
            }
            if ( id < 0 ) break;

            Point p = ray(t);
            Float density = 0.0f;
            Point4 bb;
            Point tex, q;
            Normal norm;
            Vector dpdu, dpdv;
            Vector orientation;
            Spectrum albedo;
			Float fClusterIndex = 0.f;
			int clusterIndex = 0;

            if ( m_shell.m_tetra[id].inside(m_shell.m_vtxPosition, p, bb) ) {
                const uint32_t *tmp = m_shell.m_tetra[id].idx;

                tex = m_shell.m_vtxTexcoord[tmp[0]]*bb.x +
                    m_shell.m_vtxTexcoord[tmp[1]]*bb.y +
                    m_shell.m_vtxTexcoord[tmp[2]]*bb.z +
                    m_shell.m_vtxTexcoord[tmp[3]]*bb.w;
                q = m_textureToData.transformAffine(tex);
                norm = m_shell.m_vtxNormal[tmp[0]]*bb.x +
                    m_shell.m_vtxNormal[tmp[1]]*bb.y +
                    m_shell.m_vtxNormal[tmp[2]]*bb.z +
                    m_shell.m_vtxNormal[tmp[3]]*bb.w;
                dpdu = m_shell.m_vtxTangent[tmp[0]].dpdu*bb.x +
                    m_shell.m_vtxTangent[tmp[1]].dpdu*bb.y +
                    m_shell.m_vtxTangent[tmp[2]].dpdu*bb.z +
                    m_shell.m_vtxTangent[tmp[3]].dpdu*bb.w;
                dpdv = m_shell.m_vtxTangent[tmp[0]].dpdv*bb.x +
                    m_shell.m_vtxTangent[tmp[1]].dpdv*bb.y +
                    m_shell.m_vtxTangent[tmp[2]].dpdv*bb.z +
                    m_shell.m_vtxTangent[tmp[3]].dpdv*bb.w;

                if ( m_anisotropicMedium ) {
					Frame tangFrame;
					tangFrame.s = dpdu; tangFrame.t = dpdv; tangFrame.n = norm;
					
					if (m_volume->hasOrientation()) {
						density = lookupDensity(q, ray.d, tangFrame,
							&orientation, &albedo, &fClusterIndex,
							s1, s2, pdfLobe) * m_scale;

						/*
						m_volume->lookupBundle(q, &density, &orientation, &albedo, NULL,
							NULL, NULL, &fClusterIndex);
						orientation = orientation.x*dpdu + orientation.y*dpdv + orientation.z*norm;

						if (density > 0.0f && !orientation.isZero()) {
							orientation = normalize(orientation);
							density *= m_phaseFunction->sigmaDir(dot(ray.d, orientation))*m_scale;
						}
						*/
					}
					else {
						density = lookupDensity(q, ray.d, tangFrame,
							NULL, &albedo, &fClusterIndex,
							s1, s2, pdfLobe) * m_scale;
					}
                }
                else {
					m_volume->lookupBundle(q, &density, NULL, &albedo, NULL, 
						&fClusterIndex, NULL, NULL, NULL);
                    density *= m_scale;	
                }

				if (m_useDiffAlbedoScales) {
					clusterIndex = (int)fClusterIndex;
					albedo *= m_albedoScales[clusterIndex];
				}
				else {
					clusterIndex = 0;
					albedo *= m_albedoScales[0];
				}
            }

            #if defined(HETVOL_STATISTICS)
                ++avgRayMarchingStepsSampling;
            #endif
            if (density * m_invMaxDensity > sampler->next1D()) {
                mRec.t = t;
                mRec.p = p;
                mRec.orientation = orientation;
                mRec.sigmaS = albedo * density;
				mRec.clusterIndex = clusterIndex;
				mRec.albedoScale = m_albedoScales[clusterIndex];
                mRec.sigmaA = Spectrum(density) - mRec.sigmaS;
                // XXX - what if a single channel has a 0 intensity value?
                mRec.transmittance = mRec.sigmaS.isZero() 
                    ? Spectrum(0.0f) : albedo/mRec.sigmaS;
                mRec.medium = this;
                mRec.hasExtraInfo = true;
                mRec.extra = tex;

				for (int i = 0; i < m_numLobes; i++) {
					mRec.s1[i] = s1[i];
					mRec.s2[i] = s2[i];
					mRec.pdfLobe[i] = pdfLobe[i];
					mRec.lobeScales[i] = m_lobeScales[clusterIndex][i];
				}

                success = true;
                break;
            }
        }
        mRec.medium = this;

        return success && mRec.pdfSuccess > 0;
    }


    void eval(const Ray &ray, MediumSamplingRecord &mRec) const
    {
        Log(EError, "eval(): unsupported integration method!");
    }


    bool isHomogeneous() const
    {
        return false;
    }

	const VolumeDataSourceEx *getShellmap() const {
		return m_volume.get();
	}

    std::string toString() const
    {
        std::ostringstream oss;
        oss << "HeterogeneousMediumEx3[" << endl
            << "  volume = " << indent(m_volume.toString()) << "," << endl
            << "  scale = " << m_scale << endl
            << "]";
        return oss.str();
    }


    MTS_DECLARE_CLASS()

protected:
	inline Float lookupDensity(const Point &p, const Vector &d, const Frame &tangFrame,
		Vector *_orientation = NULL, Spectrum *albedo = NULL, Float *clusterIndex = NULL,
        Spectrum *s1 = NULL, Spectrum *s2 = NULL, Float *pdfLobe = NULL) const {
		Float density;
		Vector orientation;

		Spectrum S1[MAX_SGGX_LOBES];
		Spectrum S2[MAX_SGGX_LOBES];
		Float _pdfLobe[MAX_SGGX_LOBES];
		Float weightedPdfLobe[MAX_SGGX_LOBES];

		if (_orientation) *_orientation = Vector(0.f);
		//if (pdfLobe) *pdfLobe = 0.f;
		if (clusterIndex) *clusterIndex = 0.f;

		if (m_phaseFunction->getClass()->getName() == "SGGXPhaseFunction")  {
			if (!m_volume->hasSGGXVolume()) {
				Matrix3x3 D = getPhaseFunction()->getD();
				m_volume->lookupBundle(p, &density, &orientation, albedo, NULL, 
					clusterIndex, NULL, NULL, NULL);
				orientation = orientation.x * tangFrame.s + orientation.y * tangFrame.t + orientation.z * tangFrame.n;
                
                // Update: don't need, since m_shell is in world space already
                // seems missing in original implementation
				//orientation = m_volumeToWorld(orientation);

				if (density == 0 || orientation.isZero())
					return 0.f;

				Vector w3 = normalize(orientation);
				if (_orientation) *_orientation = w3;
				Frame frame(w3);
				Matrix3x3 basis(frame.s, frame.t, w3);
				Matrix3x3 basisT;
				basis.transpose(basisT);
				Matrix3x3 S = basis * D * basisT;
				Float Sxx = S.m[0][0], Syy = S.m[1][1], Szz = S.m[2][2];
				Float Sxy = S.m[0][1], Sxz = S.m[0][2], Syz = S.m[1][2];

				S1[0][0] = Sxx; S1[0][1] = Syy; S1[0][2] = Szz;
				S2[0][0] = Sxy; S2[0][1] = Sxz; S2[0][2] = Syz;

				if (s1 && s2 && pdfLobe) {
					s1[0] = S1[0];
					s2[0] = S2[0];
					pdfLobe[0] = 1.f;
				}

				//Float sqrSum = Sxx * Sxx + Syy * Syy + Szz * Szz + Sxy * Sxy + Sxz * Sxz + Syz * Syz;
				//if (!(Sxx == 0 && Syy == 0 && Szz == 0 && Sxy == 0 && Sxz == 0 && Syz == 0))
				//if (fabsf(sqrSum) > 1e-6f)
				//	density *= m_phaseFunction->sigmaDir(d, Sxx, Syy, Szz, Sxy, Sxz, Syz);
				//else
				//	return 0.f;

				density *= m_phaseFunction->sigmaDir(d, S1, S2, pdfLobe, 1);

				return density;
			}
			else {
				bool lazy = true;

				m_volume->lookupBundle(p, &density, NULL, albedo, NULL,
					clusterIndex, S1, S2, _pdfLobe, lazy);

				if (!lazy) {
					for (int i = 0; i < m_numLobes; i++) {
						if (_pdfLobe[i] < 1e-6f) {
							S1[i] = Spectrum(0.f);
							S2[i] = Spectrum(0.f);
							continue;
						}

						Float Sxx = S1[i][0], Syy = S1[i][1], Szz = S1[i][2];
						Float Sxy = S2[i][0], Sxz = S2[i][1], Syz = S2[i][2];

						Matrix3x3 Q;
						Float eig[3];
						Matrix3x3 S;
						Vector w3, w1, w2;

						// handle orientation transform
						S = Matrix3x3(Sxx, Sxy, Sxz, Sxy, Syy, Syz, Sxz, Syz, Szz);
						S.symEig(Q, eig);
						// eig[0] < eig[1] <= eig[2]
						w3 = Vector(Q.m[0][0], Q.m[1][0], Q.m[2][0]);
						w1 = Vector(Q.m[0][1], Q.m[1][1], Q.m[2][1]);
						w2 = Vector(Q.m[0][2], Q.m[1][2], Q.m[2][2]);

						w3 = w3.x * tangFrame.s + w3.y * tangFrame.t + w3.z * tangFrame.n;

                        // Update: don't need, since m_shell is in world space already
						// seems missing in original implementation
						//w3 = m_volumeToWorld(w3);

						if (!w3.isZero()) {
							w3 = normalize(w3);

							w1 = w1.x * tangFrame.s + w1.y * tangFrame.t + w1.z * tangFrame.n;
							//w1 = m_volumeToWorld(w1);
							if (w1.isZero()) {
								Log(EInfo, "bad w1: (%.6f, %.6f, %.6f), w3: (%.6f, %.6f, %.6f)", w1.x, w1.y, w1.z,
									w3.x, w3.y, w3.z);
							}
							w1 = normalize(w1);

							w2 = w2.x * tangFrame.s + w2.y * tangFrame.t + w2.z * tangFrame.n;
							//w2 = m_volumeToWorld(w2);
							if (w2.isZero()) {
								Log(EInfo, "bad w2: (%.6f, %.6f, %.6f), w3: (%.6f, %.6f, %.6f)", w2.x, w2.y, w2.z,
									w3.x, w3.y, w3.z);
							}
							w2 = normalize(w2);

							Matrix3x3 basis(w1, w2, w3);
							Matrix3x3 D(Vector(eig[1], 0, 0), Vector(0, eig[2], 0), Vector(0, 0, eig[0]));
							Matrix3x3 basisT;
							basis.transpose(basisT);
							S = basis * D * basisT;

							Sxx = S.m[0][0]; Syy = S.m[1][1]; Szz = S.m[2][2];
							Sxy = S.m[0][1]; Sxz = S.m[0][2]; Syz = S.m[1][2];

							S1[i][0] = Sxx; S1[i][1] = Syy; S1[i][2] = Szz;
							S2[i][0] = Sxy; S2[i][1] = Sxz; S2[i][2] = Syz;
						}
						else {
							S1[i] = Spectrum(0.f);
							S2[i] = Spectrum(0.f);
						}
					}
				}
				else {
					Vector w1[MAX_SGGX_LOBES];
					Vector w2[MAX_SGGX_LOBES];
					Vector w3[MAX_SGGX_LOBES];
					Vector sigmaSqr[MAX_SGGX_LOBES];

					m_volume->lookupSGGXFrame(p, w1, w2, w3, sigmaSqr);
					
					for (int i = 0; i < m_numLobes; i++) {
						if (_pdfLobe[i] < 1e-6f) {
							S1[i] = Spectrum(0.f);
							S2[i] = Spectrum(0.f);
							continue;
						}

						w3[i] = w3[i].x * tangFrame.s + w3[i].y * tangFrame.t + w3[i].z * tangFrame.n;
						//w3[i] = m_volumeToWorld(w3[i]);

						if (!w3[i].isZero()) {
							w3[i] = normalize(w3[i]);
							w1[i] = w1[i].x * tangFrame.s + w1[i].y * tangFrame.t + w1[i].z * tangFrame.n;
							w2[i] = w2[i].x * tangFrame.s + w2[i].y * tangFrame.t + w2[i].z * tangFrame.n;
							//w1[i] = m_volumeToWorld(w1[i]);
							//w2[i] = m_volumeToWorld(w2[i]);
							w1[i] = normalize(w1[i]);
							w2[i] = normalize(w2[i]);

							Matrix3x3 basis(w1[i], w2[i], w3[i]);
							Matrix3x3 D(Vector(sigmaSqr[i][0], 0.f, 0.f), 
								Vector(0.f, sigmaSqr[i][1], 0.f), 
								Vector(0.f, 0.f, sigmaSqr[i][2]));
							Matrix3x3 basisT;
							basis.transpose(basisT);
							Matrix3x3 S = basis * D * basisT;

							S1[i][0] = S.m[0][0]; S1[i][1] = S.m[1][1]; S1[i][2] = S.m[2][2];
							S2[i][0] = S.m[0][1]; S2[i][1] = S.m[0][2]; S2[i][2] = S.m[1][2];
						}
						else {
							S1[i] = Spectrum(0.f);
							S2[i] = Spectrum(0.f);
						}
					}
				}

				if (s1 && s2 && pdfLobe) {
					for (int i = 0; i < m_numLobes; i++) {
						s1[i] = S1[i];
						s2[i] = S2[i];
						pdfLobe[i] = _pdfLobe[i];
					}
				}

				int clusterIdx = 0;
				if (*clusterIndex) clusterIdx = (int)(*clusterIndex);

				for (int i = 0; i < m_numLobes; i++) {
					weightedPdfLobe[i] = _pdfLobe[i] * m_lobeScales[clusterIdx][i];
				}
				
				density *= m_phaseFunction->sigmaDir(d, S1, S2, weightedPdfLobe, m_numLobes);

				return density;
			}
		}
		else {
			m_volume->lookupBundle(p, &density, &orientation, albedo, NULL, 
				clusterIndex, NULL, NULL, NULL);
			orientation = orientation.x * tangFrame.s + orientation.y * tangFrame.t + orientation.z * tangFrame.n;

			// seems missing in original implementation
			//orientation = m_volumeToWorld(orientation);

			if (density != 0 && !orientation.isZero()) {
				orientation = normalize(orientation);
				if (_orientation) *_orientation = orientation;
				return density*m_phaseFunction->sigmaDir(dot(d, orientation));
			}
			else
				return 0.0f;
		}
	}

protected:
    ref<VolumeDataSourceEx> m_volume;
    Float m_scale;
    bool m_anisotropicMedium;
    AABB m_volumeAABB;
    Float m_maxDensity;
    Float m_invMaxDensity;

    std::string m_shellfile;
    TetrahedronMesh m_shell;
    Transform m_worldToVolume, m_volumeToWorld, m_textureToData;

    bool m_ready;

	bool m_useDiffAlbedoScales;
	int m_numClusters;
	std::vector<Spectrum> m_albedoScales;

	int m_numLobes;
	std::vector<std::vector<Float> > m_lobeScales;
};


MTS_IMPLEMENT_CLASS_S(HeterogeneousMedium3, false, Medium)
MTS_EXPORT_PLUGIN(HeterogeneousMedium3, "Heterogeneous medium 3");
MTS_NAMESPACE_END
