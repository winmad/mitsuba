#include <mitsuba/render/volume2.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/mmap.h>
#include <mitsuba/core/math.h>


MTS_NAMESPACE_BEGIN


class InstancedVolume : public VolumeDataSourceEx {
public:
	enum EVolumeType {
		EFloat32 = 1,
		EFloat16 = 2,
		EUInt8 = 3,
		EQuantizedDirections = 4
	};


	InstancedVolume(const Properties &props) : VolumeDataSourceEx(props), m_ready(false)
    {
		m_volumeToWorld = props.getTransform("toWorld", Transform());

        m_densityFile = props.getString("densityFile");
        m_orientationFile = props.getString("orientationFile", "");

        m_indexedAlbedo = false;
        if ( props.hasProperty("albedo") )
        {
            if ( props.hasProperty("albedoFile") || props.hasProperty("indexedAlbedoFile") )
                Log(EError, "Albedo information declared more than once!");
            m_albedo = props.getSpectrum("albedo");
            m_albedoFile = "";
        }
        else if ( props.hasProperty("albedoFile") )
        {
            if ( props.hasProperty("indexedAlbedoFile") )
                Log(EError, "Albedo information declared more than once!");
            m_albedoFile = props.getString("albedoFile");
            m_albedo = Spectrum(0.0f);
        }
        else
        {
            m_albedoFile = props.getString("indexedAlbedoFile");
            m_albedo = Spectrum(0.0f);

            m_yarnColor1 = props.getSpectrum("yarnColor1");
            m_yarnColor2 = props.getSpectrum("yarnColor2");
            m_yarnSplit = props.getInteger("yarnSplit");

            m_indexedAlbedo = true;
        }

        m_blockSize = props.getVector("blockSize");
        m_divideReso.x = props.getInteger("divideX");
        m_divideReso.y = props.getInteger("divideY");

        m_reso.x = props.getInteger("tileX");
        m_reso.y = props.getInteger("tileY");
        if ( props.hasProperty("fillBlockInfo") )
        {
            if ( props.hasProperty("blockInfo") || props.hasProperty("blockFile") )
                Log(EError, "Block information declared more than once!");
            m_blockID.resize(m_reso.x*m_reso.y, props.getInteger("fillBlockInfo"));

            m_blockFile = "";
        }
        else if ( props.hasProperty("blockInfo") )
        {
            if ( props.hasProperty("blockFile") )
                Log(EError, "Block information declared more than once!");

            std::istringstream iss(props.getString("blockInfo"));
            m_blockID.resize(m_reso.x*m_reso.y);
            for ( int i = 0; i < m_reso.x*m_reso.y; ++i )
            {
                if ( !(iss >> m_blockID[i]) )
                    Log(EError, "Failed to parse the information for block %d", i);
                if ( m_blockID[i] < 0 || m_blockID[i] >= m_divideReso.x*m_divideReso.y )
                    Log(EError, "Invalid block id: %d", m_blockID[i]);
            }

            m_blockFile = "";
        }
        else
            m_blockFile = props.getString("blockFile");

        for ( int i = 0; i < 3; ++i ) m_data[i] = NULL;
	}


	InstancedVolume(Stream *stream, InstanceManager *manager) 
		: VolumeDataSourceEx(stream, manager), m_ready(false)
    {
		m_volumeToWorld = Transform(stream);

        m_densityFile = stream->readString();
        m_orientationFile = stream->readString();
        m_albedoFile = stream->readString();
        m_albedo = ( m_albedoFile == "" ? Spectrum(stream) : Spectrum(0.0f) );
        m_indexedAlbedo = ( m_albedoFile == "" ? false : stream->readBool() );
        if ( m_indexedAlbedo )
        {
            m_yarnColor1 = Spectrum(stream);
            m_yarnColor2 = Spectrum(stream);
            m_yarnSplit = stream->readInt();
        }
        m_blockSize = Vector(stream);
        m_divideReso = Vector2i(stream);

        m_reso = Vector2i(stream);
        m_blockFile = stream->readString();
        if ( m_blockFile == "" )
        {
            m_blockID.resize(m_reso.x*m_reso.y);
            stream->readIntArray(&m_blockID[0], m_blockID.size());
        }

        for ( int i = 0; i < 3; ++i ) m_data[i] = NULL;
		configure();
	}


	virtual ~InstancedVolume() {
        for ( int i = 0; i < 3; ++i )
            if ( !m_mmap[i].get() && m_data[i] ) delete m_data[i];
    }


	void serialize(Stream *stream, InstanceManager *manager) const
    {
		VolumeDataSourceEx::serialize(stream, manager);

        m_volumeToWorld.serialize(stream);

        stream->writeString(m_densityFile);
        stream->writeString(m_orientationFile);
        stream->writeString(m_albedoFile);
        if ( m_albedoFile == "" )
            m_albedo.serialize(stream);
        else
            stream->writeBool(m_indexedAlbedo);
        if ( m_indexedAlbedo )
        {

            m_yarnColor1.serialize(stream);
            m_yarnColor2.serialize(stream);
            stream->writeInt(m_yarnSplit);
        }
        m_blockSize.serialize(stream);
        m_divideReso.serialize(stream);

        m_reso.serialize(stream);
        stream->writeString(m_blockFile);
        if ( m_blockFile == "" )
            stream->writeIntArray(&m_blockID[0], m_blockID.size());
	}


	void configure()
    {
        if ( !m_ready )
        {
            m_hasOrientation = (m_orientationFile != "");

            /* Load stuff */
            loadFromFile(m_densityFile, 0, "density");
            if ( m_hasOrientation )
                loadFromFile(m_orientationFile, 1, "orientation");
            if ( m_albedoFile != "" )
                loadFromFile(m_albedoFile, 2, "albedo");

            if ( m_blockFile != "" )
            {
                ref<FileStream> fs = new FileStream(m_blockFile, FileStream::EReadOnly);
                int sz = fs->readInt();
                if ( sz != m_reso.x*m_reso.y )
                    Log(EError, "Block information size mismatch: expected %d but got %d", m_reso.x*m_reso.y, sz);

                m_blockID.resize(sz);
                fs->readIntArray(&m_blockID[0], sz);
            }

            for ( int i = 0; i < m_reso.x*m_reso.y; ++i )
                if ( m_blockID[i] < 0 || m_blockID[i] >= m_divideReso.x*m_divideReso.y )
                    Log(EError, "Invalid block id: %d", m_blockID[i]);

            /* Compute transforms */
            m_volumeAABB.min.x = -0.5f*static_cast<Float>(m_reso.x)*m_blockSize.x;
            m_volumeAABB.min.y = -0.5f*static_cast<Float>(m_reso.y)*m_blockSize.y;
            m_volumeAABB.min.z = -0.5f*m_blockSize.z;

            m_volumeAABB.max.x = -m_volumeAABB.min.x;
            m_volumeAABB.max.y = -m_volumeAABB.min.y;
            m_volumeAABB.max.z = -m_volumeAABB.min.z;

            Vector volumeExtents = m_volumeAABB.getExtents();

            m_worldToVolume = m_volumeToWorld.inverse();
            m_worldToBlock =
                Transform::scale(Vector(
                    static_cast<Float>(m_reso.x)/volumeExtents.x, 
                    static_cast<Float>(m_reso.y)/volumeExtents.y,
                    1.0f/volumeExtents.z
                ))*
                Transform::translate(Vector(
                    -m_volumeAABB.min.x, -m_volumeAABB.min.y, -m_volumeAABB.min.z
                ))*m_worldToVolume;

            m_aabb.reset();
            for ( int i = 0; i < 8; ++i )
                m_aabb.expandBy(m_volumeToWorld(m_volumeAABB.getCorner(i)));

		    /* Precompute cosine and sine lookup tables */
		    for (int i=0; i<255; i++) {
			    Float angle = (float) i * ((float) M_PI / 255.0f);
			    m_cosPhi[i] = std::cos(2.0f * angle);
			    m_sinPhi[i] = std::sin(2.0f * angle);
			    m_cosTheta[i] = std::cos(angle);
			    m_sinTheta[i] = std::sin(angle);
			    m_densityMap[i] = i/255.0f;
		    }
		    m_cosPhi[255] = m_sinPhi[255] = 0;
		    m_cosTheta[255] = m_sinTheta[255] = 0;
		    m_densityMap[255] = 1.0f;

            m_ready = true;
        }
	}


	void loadFromFile(const std::string &filename, uint32_t id, const char *dataname) {
        Log(EInfo, "Loading *%s* data from \"%s\" ..", dataname, filename.c_str());
		fs::path resolved = Thread::getThread()->getFileResolver()->resolve(filename);
        
        m_mmap[id] = new MemoryMappedFile(resolved);
        ref<MemoryStream> stream = new MemoryStream(m_mmap[id]->getData(), m_mmap[id]->getSize());
		stream->setByteOrder(Stream::ELittleEndian);

		char header[3];
		stream->read(header, 3);
		if (header[0] != 'V' || header[1] != 'O' || header[2] != 'L')
			Log(EError, "Encountered an invalid volume data file "
				"(incorrect header identifier)");
		uint8_t version;
		stream->read(&version, 1);
		if (version != 3)
			Log(EError, "Encountered an invalid volume data file "
				"(incorrect file version)");
		int type = stream->readInt();

		int xres = stream->readInt(),
			yres = stream->readInt(),
			zres = stream->readInt();
        if ( id == 0 )
		    m_dataReso = Vector3i(xres, yres, zres);
        else if ( xres != m_dataReso.x || yres != m_dataReso.y || zres != m_dataReso.z )
            Log(EError, "Specified volumes are not well-aligned!");
        m_channels[id] = stream->readInt();

		switch (type) {
			case EFloat32:
				if (m_channels[id] != 1 && m_channels[id] != 3)
					Log(EError, "Encountered an unsupported float32 volume data "
						"file (%i channels, only 1 and 3 are supported)",
						m_channels[id]);
				break;
			case EFloat16:
				Log(EError, "Error: float16 volumes are not yet supported!");
			case EUInt8:
				if (m_channels[id] != 1 && m_channels[id] != 3)
					Log(EError, "Encountered an unsupported uint8 volume data "
						"file (%i channels, only 1 and 3 are supported)", m_channels[id]);
				break;
			case EQuantizedDirections:
				if (m_channels[id] != 3)
					Log(EError, "Encountered an unsupported quantized direction "
							"volume data file (%i channels, only 3 are supported)",
							m_channels[id]);
				break;
			default:
				Log(EError, "Encountered a volume data file of unknown type (type=%i, channels=%i)!", type, m_channels);
		}
        m_volumeType[id] = static_cast<EVolumeType>(type);

		Float xmin = stream->readSingle(),
			  ymin = stream->readSingle(),
			  zmin = stream->readSingle();
		Float xmax = stream->readSingle(),
			  ymax = stream->readSingle(),
			  zmax = stream->readSingle();
        if ( id == 0 )
			m_dataAABB = AABB(Point(xmin, ymin, zmin), Point(xmax, ymax, zmax));
#if 0
        else if ( std::abs(m_dataAABB.min.x - xmin) > Epsilon || std::abs(m_dataAABB.min.y - ymin) > Epsilon || std::abs(m_dataAABB.min.z - zmin) > Epsilon ||
                  std::abs(m_dataAABB.max.x - xmax) > Epsilon || std::abs(m_dataAABB.max.y - ymax) > Epsilon || std::abs(m_dataAABB.max.z - zmax) > Epsilon )
            Log(EError, "Specified volumes are not well-aligned!");
#endif

		Log(EDebug, "Mapped \"%s\" into memory: %ix%ix%i (%i channels), %s, %s", 
			resolved.filename().c_str(), m_dataReso.x, m_dataReso.y, m_dataReso.z, m_channels[id],
			memString(m_mmap[id]->getSize()).c_str(), m_dataAABB.toString().c_str());
		m_data[id] = reinterpret_cast<uint8_t *>((reinterpret_cast<float *>(m_mmap[id]->getData())) + 12);
	}

	/**
	 * This is needed since Mitsuba might be 
	 * compiled with either single/double precision
	 */
	struct float3 {
		float value[3];

		inline float3() { }

		inline float3(float a, float b, float c) {
			value[0] = a; value[1] = b; value[2] = c;
		}

		inline float3 operator*(Float v) const {
			return float3(value[0]*v, value[1]*v, value[2]*v);
		}
		
		inline float3 operator+(const float3 &f2) const {
			return float3(value[0]+f2.value[0], value[1]+f2.value[1], value[2]+f2.value[2]);
		}

		inline Spectrum toSpectrum() const {
			Spectrum result;
			result.fromLinearRGB(value[0], value[1], value[2]);
			return result;
		}
		
		inline Vector toVector() const {
			return Vector(value[0], value[1], value[2]);
		}
	
		float operator[](int i) const {
			return value[i];
		}

		inline Matrix3x3 tensor() const {
			return Matrix3x3(
				value[0]*value[0], value[0]*value[1], value[0]*value[2],
				value[1]*value[0], value[1]*value[1], value[1]*value[2],
				value[2]*value[0], value[2]*value[1], value[2]*value[2]
			);
		}
	};


    Float lookupFloat(const Point &p) const
    {
        int idx = getCoord(p);
        if ( idx < 0 ) return 0.0;

		switch (m_volumeType[0])
        {
		case EFloat32:
            {
			    const float *floatData = (float *) m_data[0];
                return floatData[idx];
            }
		case EUInt8:
            return m_densityMap[m_data[0][idx]];
		default:
            return 0.0f;
		}
    }


    Spectrum lookupSpectrum(const Point &p) const
    {
        int idx = getCoord(p);
        if ( idx < 0 ) return Spectrum(0.0f);

        if ( m_data[2] == NULL )
            return m_albedo;
        else if ( !m_indexedAlbedo )
        {
		    switch (m_volumeType[2])
            {
			case EFloat32:
                {
				    const float3 *spectrumData = (float3 *) m_data[2];
                    return spectrumData[idx].toSpectrum();
                }
			case EUInt8:
                return float3(
                    m_densityMap[m_data[2][3*idx+0]],
                    m_densityMap[m_data[2][3*idx+1]],
                    m_densityMap[m_data[2][3*idx+2]]
                ).toSpectrum();
			default:
                return Spectrum(0.0f);
		    }
        }
        else
        {
            Spectrum ret(0.0f);

		    switch (m_volumeType[2])
            {
			case EUInt8:
                uint8_t v = m_data[2][idx];
                if ( v )
                    ret = ( v <= m_yarnSplit ? m_yarnColor1 : m_yarnColor2 );
                break;
		    }
            return ret;
        }
    }


    Vector lookupVector(const Point &p) const
    {
        Assert(m_hasOrientation);

        int idx = getCoord(p);
        if ( idx < 0 ) return Vector(0.0f);

		Vector value;
		switch (m_volumeType[1])
        {
		case EFloat32:
            {
			    const float3 *vectorData = (float3 *) m_data[1];
			    value = vectorData[idx].toVector();
			    break;
            }
		case EQuantizedDirections:
			value = lookupQuantizedDirection(idx);
			break;
		default:
            return Vector(0.0f);
		}

		if (!value.isZero())
			return normalize(m_volumeToWorld(value));
		else
			return Vector(0.0f);
    }


    void lookupBundle(const Point &p, Float *density, Vector *direction,
        Spectrum *albedo, Float *gloss) const
    {
        Assert( gloss == NULL );

        int idx = getCoord(p);
        if ( idx < 0 )
        {
            if ( density ) *density = 0.0f;
            if ( direction ) *direction = Vector(0.0f);
            if ( albedo ) *albedo = Spectrum(0.0f);
            return;
        }

        if ( density )
		    switch (m_volumeType[0])
            {
			case EFloat32:
                {
				    const float *floatData = (float *) m_data[0];
                    *density = floatData[idx];
			    }
                break;
			case EUInt8:
                *density = m_densityMap[m_data[0][idx]];
                break;
			default:
                *density = 0.0f;
		    }

        if ( direction )
        {
            Assert(m_hasOrientation);

		    Vector value;
		    switch (m_volumeType[1])
            {
		    case EFloat32:
                {
			        const float3 *vectorData = (float3 *) m_data[1];
			        value = vectorData[idx].toVector();
                }
                break;
		    case EQuantizedDirections:
			    value = lookupQuantizedDirection(idx);
			    break;
            default:
                value = Vector(0.0f);
		    }

		    if (!value.isZero())
			    *direction = normalize(m_volumeToWorld(value));
		    else
			    *direction = Vector(0.0f);
        }

        if ( albedo )
        {
            if ( m_data[2] == NULL )
                *albedo = m_albedo;
            else
            {
                *albedo = Spectrum(0.0f);

                if ( !m_indexedAlbedo )
                {
		            switch (m_volumeType[2]) {
			        case EFloat32:
                        {
				            const float3 *spectrumData = (float3 *) m_data[2];
                            *albedo = spectrumData[idx].toSpectrum();
                        }
                        break;
			        case EUInt8:
                        *albedo = float3(
                            m_densityMap[m_data[2][3*idx+0]],
                            m_densityMap[m_data[2][3*idx+1]],
                            m_densityMap[m_data[2][3*idx+2]]
                        ).toSpectrum();
                        break;
		            }
                }
                else
                {
		            switch (m_volumeType[2])
                    {
			        case EUInt8:
                        uint8_t v = m_data[2][idx];
                        if ( v )
                            *albedo = ( v <= m_yarnSplit ? m_yarnColor1 : m_yarnColor2 );
                        break;
		            }
                }
            }
        }
    }


	Float lookupFloatEx(uint32_t id, const Point &p) const
    {
        return lookupFloat(p);
	}


	Spectrum lookupSpectrumEx(uint32_t id, const Point &p) const
    {
        return lookupSpectrum(p);
	}


	Vector lookupVectorEx(uint32_t id, const Point &p) const
    {
        return lookupVector(p);
	}

    bool supportsFloatLookups() const { return true; }
	bool supportsSpectrumLookups() const { return true; }
	bool supportsVectorLookups() const { return m_hasOrientation; }
    bool supportsBundleLookups() const { return true; }
	Float getStepSize() const { return m_stepSize; }
	Float getMaximumFloatValue() const { return 1.0f; }
    Float getMaximumFloatValueEx(uint32_t id) const { return 1.0f; }

	std::string toString() const {
		std::ostringstream oss;
		oss << "InstancedVolume[" << endl
			<< "  aabb = " << m_volumeAABB.toString() << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()

protected: 
	FINLINE Vector lookupQuantizedDirection(size_t index) const {
		uint8_t theta = m_data[1][2*index], phi = m_data[1][2*index+1];
		return Vector(
			m_cosPhi[phi] * m_sinTheta[theta],
			m_sinPhi[phi] * m_sinTheta[theta],
			m_cosTheta[theta]
		);
	}

protected:
    inline int getCoord(const Point &_p) const
    {
		Point p = m_worldToBlock.transformAffine(_p);
        int bx = math::floorToInt(p.x), by = math::floorToInt(p.y);
        if ( bx < 0 || bx >= m_reso.x || by < 0 || by >= m_reso.y || p.z < 0.0f || p.z > 1.0f )
            return -1;

        int bid = m_blockID[by*m_reso.x + bx];
        int tx = bid % m_divideReso.x, ty = bid/m_divideReso.x;

        p.x = (static_cast<Float>(tx) + p.x - std::floor(p.x))/static_cast<Float>(m_divideReso.x)
            *static_cast<Float>(m_dataReso.x);
        p.y = (static_cast<Float>(ty) + p.y - std::floor(p.y))/static_cast<Float>(m_divideReso.y)
            *static_cast<Float>(m_dataReso.y);
        p.z = p.z*static_cast<Float>(m_dataReso.z);

        int ix = math::clamp(math::floorToInt(p.x), 0, m_dataReso.x - 1);
        int iy = math::clamp(math::floorToInt(p.y), 0, m_dataReso.y - 1);
        int iz = math::clamp(math::floorToInt(p.z), 0, m_dataReso.z - 1);
        return (iz*m_dataReso.y + iy)*m_dataReso.x + ix;
    }


    bool m_ready;

    std::string m_densityFile;

    std::string m_orientationFile;
    bool m_hasOrientation;

    std::string m_albedoFile;
    Spectrum m_albedo;
    bool m_indexedAlbedo;

    Spectrum m_yarnColor1, m_yarnColor2;
    int m_yarnSplit;

    Vector m_blockSize;
    Vector2i m_divideReso;
    AABB m_volumeAABB;

    std::string m_blockFile;
    Vector2i m_reso;
    std::vector<int> m_blockID;

    ref<MemoryMappedFile> m_mmap[3];
	uint8_t *m_data[3];
	EVolumeType m_volumeType[3];
    int m_channels[3];

	Vector3i m_dataReso;
	Transform m_volumeToWorld, m_worldToVolume;
    Transform m_worldToBlock;
	Float m_stepSize;
	AABB m_dataAABB;
	
	Float m_cosTheta[256], m_sinTheta[256];
	Float m_cosPhi[256], m_sinPhi[256];
	Float m_densityMap[256];
};

MTS_IMPLEMENT_CLASS_S(InstancedVolume, false, VolumeDataSourceEx);
MTS_EXPORT_PLUGIN(InstancedVolume, "Instanced volume (simple)");
MTS_NAMESPACE_END
