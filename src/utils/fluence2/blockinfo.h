#ifndef __VOLUMEINFO_H
#define __VOLUMEINFO_H


#include <mitsuba/core/aabb.h>
#include <mitsuba/core/bitmap.h>
#include "../../integrators/mft/crossimage.h"


MTS_NAMESPACE_BEGIN


class BlockInfo : public Object
{
public:
    BlockInfo(const AABB &aabb, const Vector3u &surfReso, const Vector3u &volReso);
    virtual ~BlockInfo();

    void clear();
    void addSurfacePoint(const Point &p, const Spectrum &s);
    void addVolumePoint(const Point &p, const Spectrum &s);

    void saveCrossImage(const std::string &fn) const;
    void saveVOL(const std::string &fn) const;
    //void saveVOLSlice(const std::string &fn) const;

    inline const Spectrum* getRawSurfaceData() const
    {
        return m_dataSurf;
    }
    inline const Spectrum* getRawVolumeData() const
    {
        return m_dataVol;
    }

    MTS_DECLARE_CLASS()

protected:
    AABB m_aabb, m_aabbSafe;
    Vector m_aabbSpan;
    Vector3i m_surfReso, m_volReso;

    Spectrum *m_dataSurf;
    Spectrum *m_dataVol;

    CrossImage m_cross;
};


MTS_NAMESPACE_END


#endif
