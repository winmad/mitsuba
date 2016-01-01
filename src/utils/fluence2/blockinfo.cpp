#include <mitsuba/core/fstream.h>
#include <mitsuba/core/math.h>
#include "blockinfo.h"


MTS_NAMESPACE_BEGIN


BlockInfo::BlockInfo(const AABB &aabb, const Vector3u &surfReso, const Vector3u &volReso)
    : m_aabb(aabb), m_surfReso(surfReso), m_volReso(volReso), m_cross(aabb, surfReso)
{
    m_aabbSafe = m_aabb;
    m_aabbSafe.min -= Vector(Epsilon);
    m_aabbSafe.max += Vector(Epsilon);
    m_aabbSpan = m_aabb.getExtents();

    m_dataSurf = new Spectrum[m_cross.getMaxIndex()];
    m_dataVol = new Spectrum[m_volReso.x*m_volReso.y*m_volReso.z];

    clear();
}


BlockInfo::~BlockInfo()
{
    delete[] m_dataSurf;
    delete[] m_dataVol;
}


void BlockInfo::clear()
{
    memset(m_dataSurf, 0, sizeof(Spectrum)*m_cross.getMaxIndex());
    memset(m_dataVol, 0, sizeof(Spectrum)*m_volReso.x*m_volReso.y*m_volReso.z);
}


void BlockInfo::addSurfacePoint(const Point &p, const Spectrum &s)
{
    /*
    Point2 q;
    if ( m_cross.volumeToImage(p, q) )
    {
        int qx = math::clamp(static_cast<int>(std::floor(q.x)), 0, m_dataSurfReso.x - 1);
        int qy = math::clamp(static_cast<int>(std::floor(q.y)), 0, m_dataSurfReso.y - 1);
        m_dataSurf[qy*m_dataSurfReso.x + qx] += s;
    }
    */
    uint32_t id;
    if ( m_cross.volumeToIndex(p, id) )
        m_dataSurf[id] += s;
    else
    {
        std::cout << "[Warning] bad surface point: " << p.toString() << std::endl;
        return;
    }
}


void BlockInfo::addVolumePoint(const Point &p, const Spectrum &s)
{
    if ( !m_aabbSafe.contains(p) )
    {
        std::cout << "[Warning] bad volume point: " << p.toString() << std::endl;
        return;
    }

    Point fq;
    fq.x = static_cast<float>(m_volReso.x)*(p.x - m_aabb.min.x)/m_aabbSpan.x;
    fq.y = static_cast<float>(m_volReso.y)*(p.y - m_aabb.min.y)/m_aabbSpan.y;
    fq.z = static_cast<float>(m_volReso.z)*(p.z - m_aabb.min.z)/m_aabbSpan.z;

    Point3i iq;
    iq.x = math::clamp(static_cast<int>(std::floor(fq.x)), 0, m_volReso.x - 1);
    iq.y = math::clamp(static_cast<int>(std::floor(fq.y)), 0, m_volReso.y - 1);
    iq.z = math::clamp(static_cast<int>(std::floor(fq.z)), 0, m_volReso.z - 1);

    m_dataVol[(iq.z*m_volReso.y + iq.y)*m_volReso.x + iq.x] += s;
}


void BlockInfo::saveCrossImage(const std::string &fn) const
{
    Vector2i imageSize(m_cross.getImageResolution());
    ref<Bitmap> bmp = new Bitmap(Bitmap::ERGBA, Bitmap::EFloat32, imageSize);
    Vector4* data = reinterpret_cast<Vector4*>(bmp->getFloatData());

    int index = 0;
    Point2u p;
    for ( p.y = 0; p.y < static_cast<uint32_t>(imageSize.y); ++p.y )
        for ( p.x = 0; p.x < static_cast<uint32_t>(imageSize.x); ++p.x )
        {
            Spectrum s(0.0f);
            uint32_t id;

            if ( m_cross.imageToIndex(p, id) ) s = m_dataSurf[id];
            s.toLinearRGB(data[index][0], data[index][1], data[index][2]);
            data[index++][3] = 1.0f;
            //data[index++] = Vector4(s[0], s[1], s[2], 1.0f);
        }
    ref<FileStream> stream = new FileStream(fn, FileStream::ETruncWrite);
    bmp->write(Bitmap::EOpenEXR, stream);
}


void BlockInfo::saveVOL(const std::string &fn) const
{
    ref<FileStream> stream = new FileStream(fn, FileStream::ETruncWrite);

    // Header
    stream->writeChar('V'); stream->writeChar('O'); stream->writeChar('L');
    stream->writeUChar(0x3);
    stream->writeInt(1);

    // Volume size info.
    m_volReso.serialize(stream);
    stream->writeInt(3);
    m_aabb.serialize(stream);

    // Volume data
    stream->writeFloatArray(reinterpret_cast<const Float *>(m_dataVol),
        m_volReso.x*m_volReso.y*m_volReso.z*3);
}


/*
void BlockInfo::saveVOLSlice(const std::string &fn) const
{
    ref<Bitmap> bmp = new Bitmap(Bitmap::ERGBA, Bitmap::EFloat32, Vector2i(m_volReso.x, m_volReso.z));
    Vector4* data = reinterpret_cast<Vector4*>(bmp->getFloatData());

    int index = 0;
    for ( int i = 0; i < m_volReso.z; ++i )
        for ( int j = 0; j < m_volReso.x; ++j )
        {
            const Spectrum &s = m_dataVol[((m_volReso.z - 1 - i)*m_volReso.y + m_volReso.y/2)*m_volReso.x + j];
            data[index++] = Vector4(s[0], s[1], s[2], 1.0f);
        }
    ref<FileStream> stream = new FileStream(fn, FileStream::ETruncWrite);
    bmp->write(Bitmap::EOpenEXR, stream);
}
*/


MTS_IMPLEMENT_CLASS(BlockInfo, false, Object)
MTS_NAMESPACE_END
