/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/volume2.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/fresolver.h>


MTS_NAMESPACE_BEGIN


class TiledVolumeEx : public VolumeDataSourceEx {
public:
	TiledVolumeEx(const Properties &props) : VolumeDataSourceEx(props) {
		m_volumeToWorld = props.getTransform("toWorld", Transform());
		m_stepSizeScale = props.getFloat("stepSizeScale", 1.0f);
        m_tileX = props.getInteger("tileX", 1);
        m_tileY = props.getInteger("tileY", 1);
	}

	TiledVolumeEx(Stream *stream, InstanceManager *manager) 
	: VolumeDataSourceEx(stream, manager) {
		m_volumeToWorld = Transform(stream);
		m_stepSizeScale = stream->readFloat();
        m_tileX = stream->readInt();
        m_tileY = stream->readInt();
        m_block = static_cast<VolumeDataSourceEx *>(manager->getInstance(stream));
        configure();
	}

	virtual ~TiledVolumeEx() {
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		VolumeDataSource::serialize(stream, manager);
		m_volumeToWorld.serialize(stream);
		stream->writeFloat(m_stepSizeScale);
        stream->writeInt(m_tileX);
        stream->writeInt(m_tileY);
        manager->serialize(stream, m_block.get());
	}
		
	void configure() {
        if (m_block.get() == NULL)
            Log(EError, "No embedded volume specified!");

		m_worldToVolume = m_volumeToWorld.inverse();
        m_stepSize = m_stepSizeScale*m_block->getStepSize();

        const AABB &aabb = m_block->getAABB();
        Vector extents = aabb.getExtents();
        m_blockExtents = extents;
        extents.x *= static_cast<Float>(m_tileX);
        extents.y *= static_cast<Float>(m_tileY);

        m_worldToBlock = Transform::scale(Vector(1.0f/m_blockExtents.x, 1.0f/m_blockExtents.y, 1.0f/m_blockExtents.z))*
            Transform::translate(Vector(0.5f*extents.x, 0.5f*extents.y, 0.0f))*m_worldToVolume;
        m_blockToLocal = Transform::translate(Vector(aabb.min))*Transform::scale(m_blockExtents);

        AABB aabb0(Point(-0.5f*extents.x, -0.5f*extents.y, 0.0f), Point(0.5f*extents.x, 0.5f*extents.y, extents.z));
        m_aabb.reset();
		for ( int i = 0; i < 8; ++i )
			m_aabb.expandBy(m_volumeToWorld(aabb0.getCorner(i)));

		std::ostringstream oss;
        oss << "Data AABB: " << aabb.toString() << "\nAABB: " << m_aabb.toString() << '\n';
		oss << "Step size = " << m_stepSize << " (x " << m_stepSizeScale << ")";
		Log(EDebug, oss.str().c_str());
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(VolumeDataSourceEx))) {
            Assert(m_block == NULL);
			m_block = static_cast<VolumeDataSourceEx*>(child);
		} else
			VolumeDataSource::addChild(name, child);
	}

	Float lookupFloat(const Point &_p) const {
		Point q = m_worldToBlock.transformAffine(_p);
        if ( updateLocation(q) )
            return m_block->lookupFloat(q);
        else
            return 0.0f;
	}

    Float lookupFloatEx(uint32_t id, const Point &_p) const {
		Point q = m_worldToBlock.transformAffine(_p);
        if ( updateLocation(q) )
            return m_block->lookupFloatEx(id, q);
        else
            return 0.0f;
    }

	Spectrum lookupSpectrum(const Point &_p) const {
		Point q = m_worldToBlock.transformAffine(_p);
        if ( updateLocation(q) )
            return m_block->lookupSpectrum(q);
        else
            return Spectrum(0.0f);
	}

	Spectrum lookupSpectrumEx(uint32_t id, const Point &_p) const {
		Point q = m_worldToBlock.transformAffine(_p);
        if ( updateLocation(q) )
            return m_block->lookupSpectrumEx(id, q);
        else
            return Spectrum(0.0f);
	}

	Vector lookupVector(const Point &_p) const {
		Point q = m_worldToBlock.transformAffine(_p);
        if ( updateLocation(q) ) {
            Vector ret = m_block->lookupVector(q);
            if ( !ret.isZero() ) ret = normalize(m_volumeToWorld(ret));
            return ret;
        } else
            return Vector(0.0f);
	}

	Vector lookupVectorEx(uint32_t id, const Point &_p) const {
		Point q = m_worldToBlock.transformAffine(_p);
        if ( updateLocation(q) ) {
            Vector ret = m_block->lookupVectorEx(id, q);
            if ( !ret.isZero() ) ret = normalize(m_volumeToWorld(ret));
            return ret;
        } else
            return Vector(0.0f);
	}

    void lookupBundle(const Point &_p,
        Float *density, Vector *direction, Spectrum *albedo, Float *gloss, Float *segmentation,
		Spectrum *s1, Spectrum *s2, Float *pdfLobe, bool lazy) const {
        if ( density ) *density = 0.0f;
        if ( direction ) *direction = Vector(0.0f);
        if ( albedo ) *albedo = Spectrum(0.0f);
        if ( gloss ) *gloss = 0.0f;

		for (int i = 0; i < m_block->getNumLobes(); i++) {
			if (s1) s1[i] = Spectrum(0.f);
			if (s2) s2[i] = Spectrum(0.f);
			if (pdfLobe) pdfLobe[i] = 0.f;
		}

		if (segmentation) *segmentation = 0.f;

		Point q = m_worldToBlock.transformAffine(_p);
		if ( updateLocation(q) ) {
			m_block->lookupBundle(q, density, direction, albedo, gloss, segmentation,
				s1, s2, pdfLobe, lazy);

			if (!lazy) {
				// handle orientation transform
				Matrix3x3 Q;
				Float eig[3];

				for (int i = 0; i < m_block->getNumLobes(); i++) {
					// handle orientation transform
					Matrix3x3 Q;
					Float eig[3];

					Matrix3x3 S(s1[i][0], s2[i][0], s2[i][1],
						s2[i][0], s1[i][1], s2[i][2],
						s2[i][1], s2[i][2], s1[i][2]);
					S.symEig(Q, eig);
					// eig[0] < eig[1] <= eig[2]
					Vector w3(Q.m[0][0], Q.m[1][0], Q.m[2][0]);
					w3 = m_volumeToWorld(w3);

					if (!w3.isZero()) {
						w3 = normalize(w3);

						Vector w1(Q.m[0][1], Q.m[1][1], Q.m[2][1]);
						w1 = m_volumeToWorld(w1);
						w1 = normalize(w1);
						Vector w2(Q.m[0][2], Q.m[1][2], Q.m[2][2]);
						w2 = m_volumeToWorld(w2);
						w2 = normalize(w2);

						Matrix3x3 basis(w1, w2, w3);
						Matrix3x3 D(Vector(eig[1], 0, 0), Vector(0, eig[2], 0), Vector(0, 0, eig[0]));
						Matrix3x3 basisT;
						basis.transpose(basisT);
						S = basis * D * basisT;

						s1[i][0] = S.m[0][0]; s1[i][1] = S.m[1][1]; s1[i][2] = S.m[2][2];
						s2[i][0] = S.m[0][1]; s2[i][1] = S.m[0][2]; s2[i][2] = S.m[1][2];
					}
					else {
						s1[i] = Spectrum(0.f);
						s2[i] = Spectrum(0.f);
					}
				}
			}
			else {
				for (int i = 0; i < m_block->getNumLobes(); i++) {
					s1[i] = s2[i] = Spectrum(0.f);
				}
			}
		}
    }
	
	void lookupSGGXFrame(const Point &_p,
		Vector *w1, Vector *w2, Vector *w3, Vector *sigmaSqr) const {
		Assert(w1 != NULL);
		Assert(w2 != NULL);
		Assert(w3 != NULL);
		Assert(sigmaSqr != NULL);

		Point q = m_worldToBlock.transformAffine(_p);
		if (updateLocation(q)) {
			m_block->lookupSGGXFrame(q, w1, w2, w3, sigmaSqr);

			for (int i = 0; i < m_block->getNumLobes(); i++) {
				w3[i] = m_volumeToWorld(w3[i]);

				if (!w3[i].isZero()) {
					w3[i] = normalize(w3[i]);
					w1[i] = m_volumeToWorld(w1[i]);
					w1[i] = normalize(w1[i]);
					w2[i] = m_volumeToWorld(w2[i]);
					w2[i] = normalize(w2[i]);
				}
				else {
					w1[i] = w2[i] = w3[i] = sigmaSqr[i] = Vector(0.f);
				}
			}
		}
	}

    bool supportsBundleLookups() const { return m_block->supportsBundleLookups(); }
	bool supportsFloatLookups() const { return m_block->supportsFloatLookups(); }
    bool supportsSpectrumLookups() const { return m_block->supportsSpectrumLookups(); }
	bool supportsVectorLookups() const { return m_block->supportsVectorLookups(); }
	Float getStepSize() const { return m_stepSize; }
    Float getMaximumFloatValue() const { return m_block->getMaximumFloatValue(); }
    Float getMaximumFloatValueEx(uint32_t id) const { return m_block->getMaximumFloatValueEx(id); }

	bool hasOrientation() const {
		return m_block->hasOrientation();
	}

	bool hasSGGXVolume() const {
		return m_block->hasSGGXVolume();
	}

	int getNumLobes() const {
		return m_block->getNumLobes();
	}

	MTS_DECLARE_CLASS()

protected:
    inline bool updateLocation(Point &q) const {
        if ( q.x > -Epsilon && q.x < static_cast<Float>(m_tileX) + Epsilon && 
             q.y > -Epsilon && q.y < static_cast<Float>(m_tileY) + Epsilon &&
             q.z > -Epsilon && q.z < 1.0f + Epsilon ) {
            q.x -= std::floor(q.x); q.y -= std::floor(q.y);
            q = m_blockToLocal.transformAffine(q);
            return true;
        } else
            return false;
    }

    ref<VolumeDataSourceEx> m_block;
	Transform m_worldToVolume, m_volumeToWorld;
    Transform m_worldToBlock, m_blockToLocal;
	Float m_stepSize, m_stepSizeScale;

    Vector m_blockExtents;
    int m_tileX, m_tileY;
};

MTS_IMPLEMENT_CLASS_S(TiledVolumeEx, false, VolumeDataSourceEx);
MTS_EXPORT_PLUGIN(TiledVolumeEx, "Tiled volume data source 2");
MTS_NAMESPACE_END
