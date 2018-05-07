#pragma once
#if !defined(__SPHERICAL_DISTRIBUTION_H_)
#define  __SPHERICAL_DISTRIBUTION_H_

#include <mitsuba/core/sched.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <boost/filesystem/path.hpp>

#define USE_SQUARE_CONCENTRIC

MTS_NAMESPACE_BEGIN

class MTS_EXPORT_RENDER SphericalDistribution : public WorkResult {
public:
	SphericalDistribution(int size);
	void clear();
	void put(const SphericalDistribution *dist, bool putBitmap=true);
	void put(const Vector &dir, const Spectrum &value, double weight, double normFactor, bool putBitmap=true);
	void scale(double scale);
	void saveExr(fs::path filename);

	void load(Stream *stream);
	void save(Stream *stream) const;
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	virtual ~SphericalDistribution() {}
public:
	int m_size;
	Vector3d m_totValue;
	int m_totValidParticles;
	double m_totWeight;

	ref<Bitmap> m_values;
	
// 	Vector3d m_moments[3];
// 	Float cY[3];
};

class MTS_EXPORT_RENDER MultiLobeDistribution : public WorkResult {
public:
	MultiLobeDistribution(int numLobes, int size);
	SphericalDistribution *getLobe(int lobeIdx);
	const SphericalDistribution *getLobe(int lobeIdx) const;
	void clear();
	void put(const MultiLobeDistribution *dist);

	void load(Stream *stream);
	void save(Stream *stream) const;
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	virtual ~MultiLobeDistribution() {}
public:
	int m_numLobes;
	std::vector<ref<SphericalDistribution> > m_lobes;
};

MTS_NAMESPACE_END

#endif
