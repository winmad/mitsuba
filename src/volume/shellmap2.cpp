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

#define TETRAHEDRON_MESH_NO_CACHE
#include "tetra.h"


MTS_NAMESPACE_BEGIN


class ShellMappedDataSourceEx : public VolumeDataSourceEx {
public:
	ShellMappedDataSourceEx(const Properties &props) : VolumeDataSourceEx(props) {
		m_volumeToWorld = props.getTransform("toWorld", Transform());
		m_shellfile = props.getString("shellfile");
		fs::path resolved = Thread::getThread()->getFileResolver()->resolve(m_shellfile);
		if ( !m_shell.load(resolved.string().c_str()) )
			Log(EError, "Failed to load the shell file!");
        else
            Log(EInfo, "Shell mesh loaded: %u tetrahedra, tree depth: %u",
                m_shell.getTetrahedronCount(), m_shell.getTreeDepth());
		m_stepSizeScale = props.getFloat("stepSizeScale", 1.0f);
	}

	ShellMappedDataSourceEx(Stream *stream, InstanceManager *manager) 
	: VolumeDataSourceEx(stream, manager) {
		m_volumeToWorld = Transform(stream);
		m_shellfile = stream->readString();
		if ( !m_shell.load(m_shellfile.c_str()) )
			Log(EError, "Failed to load the shell file!");
        else
            Log(EInfo, "Shell mesh loaded: %u tetrahedra, tree depth: %u",
                m_shell.getTetrahedronCount(), m_shell.getTreeDepth());
		m_stepSizeScale = stream->readFloat();
        m_block = static_cast<VolumeDataSourceEx *>(manager->getInstance(stream));
        configure();
	}

	virtual ~ShellMappedDataSourceEx() {
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		VolumeDataSource::serialize(stream, manager);
		m_volumeToWorld.serialize(stream);
		stream->writeString(m_shellfile);
		stream->writeFloat(m_stepSizeScale);
        manager->serialize(stream, m_block.get());
	}
		
	void configure() {
        if (m_block.get() == NULL)
            Log(EError, "No embedded volume specified!");

		m_worldToVolume = m_volumeToWorld.inverse();
        const AABB &aabb = m_block->getAABB();
        m_textureToData = Transform::translate(Vector(aabb.min))*Transform::scale(aabb.getExtents());

        m_stepSize = m_stepSizeScale*m_block->getStepSize();
        m_aabb.reset();
		for ( int i = 0; i < 8; ++i )
			m_aabb.expandBy(m_volumeToWorld(m_shell.getAABB().getCorner(i)));

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

private:
    void clampPoint(Point &p) const {
        Assert( p.z > -Epsilon && 1.0f - p.z > -Epsilon );
        p.x -= std::floor(p.x);
        p.y -= std::floor(p.y);
    }
    
public:
	Float lookupFloat(const Point &_p) const {
		Point q = m_worldToVolume.transformAffine(_p);
		Point p;
		if ( !m_shell.lookupPoint(q, p) ) return 0.0f;
        clampPoint(p);
        return m_block->lookupFloat(m_textureToData.transformAffine(p));
	}

	Float lookupFloatEx(uint32_t id, const Point &_p) const {
		Point q = m_worldToVolume.transformAffine(_p);
		Point p;
		if ( !m_shell.lookupPoint(q, p) ) return 0.0f;
        clampPoint(p);
        return m_block->lookupFloatEx(id, m_textureToData.transformAffine(p));
	}
    
	Spectrum lookupSpectrum(const Point &_p) const {
		Point q = m_worldToVolume.transformAffine(_p);
		Point p;
		if ( !m_shell.lookupPoint(q, p) ) return Spectrum(0.0f);
        clampPoint(p);
        return m_block->lookupSpectrum(m_textureToData.transformAffine(p));
	}

	Spectrum lookupSpectrumEx(uint32_t _id, const Point &_p) const {
		Point q = m_worldToVolume.transformAffine(_p);
		Point p;
		if ( !m_shell.lookupPoint(q, p) ) return Spectrum(0.0f);
        clampPoint(p);
        return m_block->lookupSpectrumEx(_id, m_textureToData.transformAffine(p));
	}

	Vector lookupVector(const Point &_p) const {
		Point q = m_worldToVolume.transformAffine(_p);
		Point p;
		Vector norm;
		TangentSpace tang;
		if ( !m_shell.lookupPoint(q, p, norm, tang) ) return Vector(0.0f);
        clampPoint(p);

        Vector tmp = m_block->lookupVector(m_textureToData.transformAffine(p));
        Vector ret = m_volumeToWorld(tmp.x*tang.dpdu + tmp.y*tang.dpdv + tmp.z*norm);
		if ( !ret.isZero() ) ret = normalize(ret);
		return ret;
	}

	Vector lookupVectorEx(uint32_t _id, const Point &_p) const {
		Point q = m_worldToVolume.transformAffine(_p);
		Point p;
		Vector norm;
		TangentSpace tang;
		if ( !m_shell.lookupPoint(q, p, norm, tang) ) return Vector(0.0f);
        clampPoint(p);

        Vector tmp = m_block->lookupVectorEx(_id, m_textureToData.transformAffine(p));
        Vector ret = m_volumeToWorld(tmp.x*tang.dpdu + tmp.y*tang.dpdv + tmp.z*norm);
		if ( !ret.isZero() ) ret = normalize(ret);
		return ret;
	}

	// doesn't implement lazy SGGX evaluation
	void lookupBundle(const Point &_p, Float *density, Vector *direction,
		Spectrum *albedo, Float *gloss, Float *segmentation,
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

		Point q = m_worldToVolume.transformAffine(_p);
		Point p;

        if ( direction ) {
		    Vector norm;
		    TangentSpace tang;
            if ( m_shell.lookupPoint(q, p, norm, tang) ) {
                clampPoint(p);
                m_block->lookupBundle(m_textureToData.transformAffine(p), density, direction, albedo, gloss,
					segmentation, s1, s2, pdfLobe);
                *direction = m_volumeToWorld(direction->x*tang.dpdu + direction->y*tang.dpdv + direction->z*norm);
		        if ( !direction->isZero() ) *direction = normalize(*direction);
            }
        }
		else if (s1) {
			Vector norm;
			TangentSpace tang;
			if (m_shell.lookupPoint(q, p, norm, tang)) {
				clampPoint(p);
				m_block->lookupBundle(m_textureToData.transformAffine(p), density, direction, albedo, gloss,
					segmentation, s1, s2, pdfLobe, lazy);

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
						w3 = m_volumeToWorld(w3.x * tang.dpdu + w3.y * tang.dpdv + w3.z * norm);

						if (!w3.isZero()) {
							w3 = normalize(w3);

							Vector w1(Q.m[0][1], Q.m[1][1], Q.m[2][1]);
							w1 = m_volumeToWorld(w1.x * tang.dpdu + w1.y * tang.dpdv + w1.z * norm);
							w1 = normalize(w1);
							Vector w2(Q.m[0][2], Q.m[1][2], Q.m[2][2]);
							w2 = m_volumeToWorld(w2.x * tang.dpdu + w2.y * tang.dpdv + w2.z * norm);
							w2 = normalize(w2);

							//Frame frame(w3);
							//Matrix3x3 basis(frame.s, frame.t, w3);

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
		else {
            if ( m_shell.lookupPoint(q, p) ) {
                clampPoint(p);
                m_block->lookupBundle(m_textureToData.transformAffine(p), density, NULL, albedo, gloss,
					segmentation, s1, s2, pdfLobe);
            }
        }
    }


	void lookupSGGXFrame(const Point &_p,
		Vector *w1, Vector *w2, Vector *w3, Vector *sigmaSqr) const {
		Assert(w1 != NULL);
		Assert(w2 != NULL);
		Assert(w3 != NULL);
		Assert(sigmaSqr != NULL);

		Point q = m_worldToVolume.transformAffine(_p);
		Point p;

		Vector norm;
		TangentSpace tang;
		if (m_shell.lookupPoint(q, p, norm, tang)) {
			clampPoint(p);
			m_block->lookupSGGXFrame(m_textureToData.transformAffine(p),
				w1, w2, w3, sigmaSqr);

			for (int i = 0; i < m_block->getNumLobes(); i++) {
				w3[i] = m_volumeToWorld(w3[i].x * tang.dpdu + w3[i].y * tang.dpdv + w3[i].z * norm);

				if (!w3[i].isZero()) {
					w3[i] = normalize(w3[i]);
					w1[i] = m_volumeToWorld(w1[i].x * tang.dpdu + w1[i].y * tang.dpdv + w1[i].z * norm);
					w1[i] = normalize(w1[i]);
					w2[i] = m_volumeToWorld(w2[i].x * tang.dpdu + w2[i].y * tang.dpdv + w2[i].z * norm);
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
	std::string m_shellfile;
    TetrahedronMesh m_shell;
    ref<VolumeDataSourceEx> m_block;
	Transform m_worldToVolume, m_volumeToWorld;
    Transform m_textureToData;
	Float m_stepSize, m_stepSizeScale;
};

MTS_IMPLEMENT_CLASS_S(ShellMappedDataSourceEx, false, VolumeDataSourceEx);
MTS_EXPORT_PLUGIN(ShellMappedDataSourceEx, "Shell-mapped data source 2");
MTS_NAMESPACE_END
