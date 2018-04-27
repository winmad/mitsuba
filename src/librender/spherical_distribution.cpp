#include <mitsuba/render/spherical_distribution.h>

MTS_NAMESPACE_BEGIN

// WorkResult
SphericalDistribution::SphericalDistribution(int size) : m_size(size) {
	m_values = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat64, Vector2i(m_size));
// 	cY[0] = 0.5 * std::sqrt(INV_PI);
// 	cY[1] = 0.5 * std::sqrt(3.0 * INV_PI);
// 	cY[2] = 0.25 * std::sqrt(5.0 * INV_PI);
}

void SphericalDistribution::clear() {
	m_totValue = Vector3d(0.0);
	m_totValidParticles = 0;
	m_totWeight = 0.0;
	double *data = m_values->getFloat64Data();
	for (int i = 0; i < m_size * m_size * SPECTRUM_SAMPLES; i++)
		*data++ = 0.0;

// 	for (int i = 0; i < 3; i++)
// 		m_moments[i] = Vector3d(0.0f);
}

void SphericalDistribution::put(const SphericalDistribution *dist, bool putBitmap) {
	Assert(m_size == dist->m_size);
	m_totValue += dist->m_totValue;
	m_totValidParticles += dist->m_totValidParticles;
	m_totWeight += dist->m_totWeight;

// 	for (int i = 0; i < 3; i++)
// 		m_moments[i] += dist->m_moments[i];

	if (putBitmap) {
		double *target = m_values->getFloat64Data();
		const double *source = dist->m_values->getFloat64Data();
		for (int i = 0; i < m_size * m_size * SPECTRUM_SAMPLES; i++)
			*target++ += *source++;
	}
}

void SphericalDistribution::put(const Vector &dir, const Spectrum &value, double weight, double normFactor, bool putBitmap) {
	for (int k = 0; k < 3; k++) {
		m_totValue[k] += (double)value[k] * weight;
	}
	m_totValidParticles++;
	m_totWeight += weight;

// 	for (int i = 0; i < 3; i++) {
// 		Float w = cY[i];
// 		if (i == 1) {
// 			w *= dir.z;
// 		}
// 		else if (i == 2) {
// 			w *= (3.0 * dir.z * dir.z - 1.0);
// 		}
// 
// 		for (int c = 0; c < 3; c++) {
// 			m_moments[i][c] += (double)value[c] * w;
// 		}
// 	}

	if (putBitmap) {
		int c = math::clamp(math::floorToInt((dir.x + 1.0) * 0.5 * m_size), 0, m_size - 1);
		int r = math::clamp(math::floorToInt((dir.y + 1.0) * 0.5 * m_size), 0, m_size - 1);

		double *data = m_values->getFloat64Data();
		int idx = (m_size - r - 1) * m_size + c;
		for (int k = 0; k < 3; k++) {
			data[3 * idx + k] += (double)value[k] * weight * (dir.z * normFactor);
		}
	}
}

void SphericalDistribution::scale(double scale) {
	double *data = m_values->getFloat64Data();
	for (int i = 0; i < m_size * m_size * SPECTRUM_SAMPLES; i++)
		*data++ *= scale;
}

void SphericalDistribution::saveExr(fs::path filename) {
	filename.replace_extension(".exr");
	ref<Bitmap> bitmap;
	bitmap = m_values->convert(Bitmap::ERGB, Bitmap::EFloat32);
	ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
	bitmap->write(Bitmap::EOpenEXR, stream);
}

void SphericalDistribution::load(Stream *stream) {
	m_size = stream->readInt();
	m_totValue[0] = stream->readDouble();
	m_totValue[1] = stream->readDouble();
	m_totValue[2] = stream->readDouble();
	m_totValidParticles = stream->readInt();
	m_totWeight = stream->readDouble();
	stream->readDoubleArray(m_values->getFloat64Data(), m_size * m_size * SPECTRUM_SAMPLES);
// 	for (int i = 0; i < 3; i++)
// 		for (int c = 0; c < 3; c++)
// 			m_moments[i][c] = stream->readDouble();
}

void SphericalDistribution::save(Stream *stream) const {
	stream->writeInt(m_size);
	stream->writeDouble(m_totValue[0]);
	stream->writeDouble(m_totValue[1]);
	stream->writeDouble(m_totValue[2]);
	stream->writeInt(m_totValidParticles);
	stream->writeDouble(m_totWeight);
	stream->writeDoubleArray(m_values->getFloat64Data(), m_size * m_size * SPECTRUM_SAMPLES);
// 	for (int i = 0; i < 3; i++)
// 		for (int c = 0; c < 3; c++)
// 			stream->writeDouble(m_moments[i][c]);
}

std::string SphericalDistribution::toString() const {
	std::ostringstream oss;
	oss << "SphericalDistribution[" << endl
		<< "  size = " << m_size << endl
		<< "]";
	return oss.str();
}

// WorkResult
MultiLobeDistribution::MultiLobeDistribution(int numLobes, int size) : m_numLobes(numLobes) {
	m_lobes.resize(numLobes);
	for (int i = 0; i < numLobes; i++) {
		m_lobes[i] = new SphericalDistribution(size);
	}
}

SphericalDistribution* MultiLobeDistribution::getLobe(int lobeIdx) {
	return m_lobes[lobeIdx].get();
}

const SphericalDistribution* MultiLobeDistribution::getLobe(int lobeIdx) const {
	return m_lobes[lobeIdx].get();
}

void MultiLobeDistribution::clear() {
	for (int i = 0; i < m_numLobes; i++) {
		m_lobes[i]->clear();
	}
}

void MultiLobeDistribution::put(const MultiLobeDistribution *dist) {
	for (int i = 0; i < m_numLobes; i++) {
		m_lobes[i]->put(dist->getLobe(i));
	}
}

void MultiLobeDistribution::load(Stream *stream) {
	m_numLobes = stream->readInt();
	m_lobes.resize(m_numLobes);
	for (int i = 0; i < m_numLobes; i++) {
		m_lobes[i]->load(stream);
	}
}

void MultiLobeDistribution::save(Stream *stream) const {
	stream->writeInt(m_numLobes);
	for (int i = 0; i < m_numLobes; i++) {
		m_lobes[i]->save(stream);
	}
}

std::string MultiLobeDistribution::toString() const {
	std::ostringstream oss;
	oss << "MultiLobeDistribution[" << endl
		<< "  numLobes = " << m_numLobes << endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(SphericalDistribution, false, WorkResult)
MTS_IMPLEMENT_CLASS(MultiLobeDistribution, false, WorkResult)
MTS_NAMESPACE_END
