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

#if !defined(__VOLUME_EX_H)
#define __VOLUME_EX_H

#include <mitsuba/render/volume.h>

//#define USE_STOC_EVAL // using stochastic evaluation of SGGX

MTS_NAMESPACE_BEGIN

class MTS_EXPORT_RENDER VolumeDataSourceEx : public VolumeDataSource {
public:
    virtual Float lookupFloatEx(uint32_t id, const Point &p) const;
    virtual Spectrum lookupSpectrumEx(uint32_t id, const Point &p) const;
    virtual Vector lookupVectorEx(uint32_t id, const Point &p) const;

    virtual bool supportsBundleLookups() const;
    virtual void lookupBundle(const Point &p,
        Float *density, Vector *direction, Spectrum *albedo, Float *gloss, Float *segmentation, 
		Spectrum *s1, Spectrum *s2, Float *pdfLobe, bool lazy = false) const;

    virtual Float getMaximumFloatValueEx(uint32_t id) const = 0;

	virtual bool hasOrientation() const;
	virtual bool hasSGGXVolume() const;

	virtual int getNumLobes() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~VolumeDataSourceEx();

	/// Protected constructor
	VolumeDataSourceEx(const Properties &props);

	/// Unserialize from a binary data stream
	VolumeDataSourceEx(Stream *stream, InstanceManager *manager); 
};

MTS_NAMESPACE_END

#endif /* __VOLUME_EX_H */
