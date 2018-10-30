#include <mitsuba/render/vectorblock.h>

MTS_NAMESPACE_BEGIN

VectorBlock::VectorBlock(int dimensions) {
	m_dimensions = dimensions;
	m_data.resize(m_dimensions);
	clear();
}

void VectorBlock::load(Stream *stream) {
	m_offset = Point2i(stream);
	m_dimensions = stream->readInt();
	m_data.resize(m_dimensions);
	stream->readFloatArray(m_data.data(), m_dimensions);
	m_weight = stream->readFloat();
}

void VectorBlock::save(Stream *stream) const {
	m_offset.serialize(stream);
	stream->writeInt(m_dimensions);
	stream->writeFloatArray(m_data.data(), m_data.size());
	stream->writeFloat(m_weight);
}

std::string VectorBlock::toString() const {
	std::ostringstream oss;
	oss << "VectorBlock[" << endl
		<< "  dimensions = " << m_dimensions << endl
		<< "  weight = " << m_weight << endl
		<< "  data = [";
	for (auto i = m_data.begin(); i != m_data.end(); ++i)
		oss << *i << ", ";
	oss << "]" << endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(VectorBlock, false, WorkResult)
MTS_NAMESPACE_END
