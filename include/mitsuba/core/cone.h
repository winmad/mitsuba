/*
	Added by Lifan Wu
	Oct 18, 2015
*/

#pragma once
#if !defined(__MITSUBA_CORE_CONE_H_)
#define __MITSUBA_CORE_CONE_H_

#include <mitsuba/core/aabb.h>
#include <mitsuba/core/math.h>

MTS_NAMESPACE_BEGIN

struct Cone {
	inline Cone() {}
	inline Cone(const Vector &v) : dirBox((Point) v) {}
	inline Cone(const Cone& c) : dirBox(c.dirBox) {}

	inline void reset() {
		dirBox.reset();
	}

	inline Vector getAxis() const {
		Vector result(dirBox.getCenter());
		return normalize(result);
	}

	inline Float getAngleCos() const {
		Vector center = Vector(dirBox.getCenter());
		Float r2 = (dirBox.max - dirBox.min).lengthSquared() * 0.25f;
		Float d2 = center.lengthSquared();
		if (d2 == 0) return 0.f;
		Float d = sqrt(d2);
		Float result = (d2 - r2 + 1.f) / (d * 2.f);
		return std::min(std::max(result, 0.f), 1.f);
	}

	inline Float getAngleSin() const {
		Float angleCos = getAngleCos();
		return sqrt(1.f - angleCos * angleCos);
	}

	inline void build(const std::vector<Vector> &v) {
		for (size_t i = 0; i < v.size(); i++) {
			dirBox.expandBy((Point) v[i]);
		}
	}

	inline bool isValid() const {
		Vector center = Vector(dirBox.getCenter());
		Float r2 = (dirBox.max - dirBox.min).lengthSquared() * 0.25f;
		Float d2 = center.lengthSquared();
		if (d2 == 0) return false;
		Float d = sqrt(d2);
		Float angleCos = (d2 - r2 + 1.f) / (d * 2.f);
		return (angleCos >= 0.f && angleCos <= 1.f);
	}

	inline bool contains(const Vector &v) const {
		return dirBox.contains((Point) v);
	}

	inline bool overlaps(const Cone &c) const {
		for (int i=0; i<3; ++i)
			if (dirBox.max[i] < c.dirBox.min[i] || dirBox.min[i] > c.dirBox.max[i])
				return false;
		return true;
	}

	inline void expandBy(const Vector &v) {
		dirBox.expandBy((Point) v);
	}

	inline void expandBy(const Cone &c) {
		dirBox.expandBy(c.dirBox);
	}

	inline bool operator==(const Cone &c) const {
		return dirBox == c.dirBox;
	}

	inline bool operator!=(const Cone &c) const {
		return dirBox != c.dirBox;
	}

	AABB dirBox;
};

MTS_NAMESPACE_END

#endif