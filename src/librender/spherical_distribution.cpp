#include <mitsuba/render/spherical_distribution.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

// WorkResult
SphericalDistribution::SphericalDistribution(int size, int useFullSphere) : m_size(size), m_useFullSphere(useFullSphere) {
	if (m_useFullSphere)
		m_imgSize = Vector2i(m_size, m_size * 2);
	else
		m_imgSize = Vector2i(m_size);
	
	m_values = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat64, m_imgSize);

// 	cY[0] = 0.5 * std::sqrt(INV_PI);
// 	cY[1] = 0.5 * std::sqrt(3.0 * INV_PI);
// 	cY[2] = 0.25 * std::sqrt(5.0 * INV_PI);
}

void SphericalDistribution::clear() {
	m_totValue = Vector3d(0.0);
	m_totValidParticles = 0;
	m_totWeight = 0.0;
	double *data = m_values->getFloat64Data();

	for (int i = 0; i < m_imgSize.x * m_imgSize.y * SPECTRUM_SAMPLES; i++)
		*data++ = 0.0;

// 	for (int i = 0; i < 3; i++)
// 		m_moments[i] = Vector3d(0.0f);
}

void SphericalDistribution::put(const SphericalDistribution *dist, bool putBitmap) {
	Assert(m_imgSize == dist->m_imgSize);
	m_totValue += dist->m_totValue;
	m_totValidParticles += dist->m_totValidParticles;
	m_totWeight += dist->m_totWeight;

// 	for (int i = 0; i < 3; i++)
// 		m_moments[i] += dist->m_moments[i];

	if (putBitmap) {
		double *target = m_values->getFloat64Data();
		const double *source = dist->m_values->getFloat64Data();
		for (int i = 0; i < m_imgSize.x * m_imgSize.y * SPECTRUM_SAMPLES; i++)
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
#if defined(USE_SQUARE_CONCENTRIC)
		int offset_r = 0;
		if (dir.z <= 0) {
			if (m_useFullSphere)
				offset_r = m_size;
			else
				return;
		}

		double *data = m_values->getFloat64Data();
		Point2 square = warp::uniformHemisphereToSquareConcentric(dir);

		/*
		int c = math::clamp(math::floorToInt(square.x * m_size), 0, m_size - 1);
		int r = math::clamp(math::floorToInt(square.y * m_size), 0, m_size - 1);
		
		for (int k = 0; k < 3; k++) {
			int idx = r * m_size + c;
			data[3 * idx + k] += value[k] * weight * normFactor;
		}
		*/

		int numCells = m_size - 1;
		int c = math::clamp(math::floorToInt(square.x * numCells), 0, numCells - 1);
		int r = math::clamp(math::floorToInt(square.y * numCells), 0, numCells - 1);

		Float u = square.x * numCells - c;
		Float v = square.y * numCells - r;
		
		Spectrum tmp = value * weight * normFactor;
		// small kernel
		for (int dr = 0; dr < 2; dr++) {
			double wv = std::abs(1.0 - dr - v);
			for (int dc = 0; dc < 2; dc++) {
				double wKernel = wv * std::abs(1.0 - dc - u);
				int idx = (r + dr + offset_r) * m_size + (c + dc);
				for (int k = 0; k < 3; k++) {
					data[3 * idx + k] += tmp[k] * wKernel;
				}
			}
		}
		
#else
		int c = math::clamp(math::floorToInt((dir.x + 1.0) * 0.5 * m_size), 0, m_size - 1);
		int r = math::clamp(math::floorToInt((dir.y + 1.0) * 0.5 * m_size), 0, m_size - 1);

		double *data = m_values->getFloat64Data();
		int idx = (m_size - r - 1) * m_size + c;
		for (int k = 0; k < 3; k++) {
			data[3 * idx + k] += (double)value[k] * weight * (dir.z * normFactor);
		}
#endif
	}
}

void SphericalDistribution::scale(double scale) {
	double *data = m_values->getFloat64Data();
	for (int i = 0; i < m_imgSize.x * m_imgSize.y * SPECTRUM_SAMPLES; i++)
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
	m_useFullSphere = stream->readInt();
	m_imgSize.x = stream->readInt();
	m_imgSize.y = stream->readInt();
	m_totValue[0] = stream->readDouble();
	m_totValue[1] = stream->readDouble();
	m_totValue[2] = stream->readDouble();
	m_totValidParticles = stream->readInt();
	m_totWeight = stream->readDouble();
	stream->readDoubleArray(m_values->getFloat64Data(), m_imgSize.x * m_imgSize.y * SPECTRUM_SAMPLES);
// 	for (int i = 0; i < 3; i++)
// 		for (int c = 0; c < 3; c++)
// 			m_moments[i][c] = stream->readDouble();
}

void SphericalDistribution::save(Stream *stream) const {
	stream->writeInt(m_size);
	stream->writeInt(m_useFullSphere);
	stream->writeInt(m_imgSize.x);
	stream->writeInt(m_imgSize.y);
	stream->writeDouble(m_totValue[0]);
	stream->writeDouble(m_totValue[1]);
	stream->writeDouble(m_totValue[2]);
	stream->writeInt(m_totValidParticles);
	stream->writeDouble(m_totWeight);
	stream->writeDoubleArray(m_values->getFloat64Data(), m_imgSize.x * m_imgSize.y * SPECTRUM_SAMPLES);
// 	for (int i = 0; i < 3; i++)
// 		for (int c = 0; c < 3; c++)
// 			stream->writeDouble(m_moments[i][c]);
}

std::string SphericalDistribution::toString() const {
	std::ostringstream oss;
	oss << "SphericalDistribution[" << endl
		<< "  size = " << m_size << endl
		<< "  useFullSphere = " << m_useFullSphere << endl
		<< "]";
	return oss.str();
}

// WorkResult
MultiLobeDistribution::MultiLobeDistribution(int numLobes, int size, int useFullSphere) : m_numLobes(numLobes) {
	m_lobes.resize(numLobes);
	for (int i = 0; i < numLobes; i++) {
		m_lobes[i] = new SphericalDistribution(size, useFullSphere);
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
