#pragma once
#if !defined(__MITSUBA_RENDER_VECTORBLOCK_H_)
#define __MITSUBA_RENDER_VECTORBLOCK_H_

#include <mitsuba/core/sched.h>

MTS_NAMESPACE_BEGIN

class MTS_EXPORT_RENDER VectorBlock : public WorkResult {
public:
    VectorBlock(int dimensions);

    /// Set the current block offset
	inline void setOffset(const Point2i &offset) { m_offset = offset; }
    
    /// Return the current block offset
    inline const Point2i &getOffset() const { return m_offset; }

    inline int getDimensions() const { return m_dimensions; }

    inline void setDimensions(int dimensions) {
        m_dimensions = dimensions;
        m_data.resize(m_dimensions);
        clear();
    }

    inline const std::vector<Float> &getData() const { return m_data; }

    inline double getDataAt(int idx) const { return m_data[idx]; }

    inline void setData(const std::vector<Float> &data) { m_data = data; }

    inline double getWeight() const { return m_weight; }

    inline void setWeight(double weight) { m_weight = weight; }

	inline void addWeight(double weight) { m_weight += weight; }

    /*
    inline uint8_t *getUInt8Data() {
        double *data = (double *) m_dataBytes;
        for (int i = 0; i < m_dimensions; i++) {
            data[i] = m_data[i];
        }
        return m_dataBytes;
    }
    */

	/// Clear everything to zero
	inline void clear() { 
        memset(&m_data[0], 0, m_dimensions * sizeof(Float));
        m_weight = 0.0f;
    }

    inline void put(const VectorBlock *block) {
        for (int i = 0; i < m_dimensions; i++)
            m_data[i] += block->m_data[i];
        m_weight += block->m_weight;
    }

    inline bool put(const Float *value, const Float w) {
        // check if all sample values are valid
        for (int i = 0; i < m_dimensions; i++) {
            if (EXPECT_NOT_TAKEN(!std::isfinite(value[i]))) {
                std::ostringstream oss;
                oss << "Invalid sample value : [";
                for (int j = 0; j < m_dimensions; j++) {
                    oss << value[j];
                    if (j + 1 < m_dimensions)
                        oss << ", ";
                }
                oss << "]";
                Log(EWarn, "%s", oss.str().c_str());
                return false;
            }
        }

        for (int i = 0; i < m_dimensions; i++) {
            m_data[i] += value[i];
        }
        m_weight += w;
        return true;
    }

    inline bool put(const std::vector<int>& indices, const std::vector<Float>& values, const Float w) {
        // check if all sample values are valid
        for (int i = 0; i < values.size(); i++) {
			Float value = values[i];
            if (EXPECT_NOT_TAKEN(!std::isfinite(value))) {
                std::ostringstream oss;
                oss << "Invalid sample value : [" ;
                for (int j = 0, jend = values.size(); j < jend; j++) {
                    oss << indices[j] << ": " << values[j];
                    if (j + 1 < jend)
                       oss << ", ";
                }
                oss << "]";
                Log(EWarn, "%s", oss.str().c_str());
                return false;
            }
        }

        for (int i = 0, iend = values.size(); i < iend; i++) {
            m_data[indices[i]] += values[i];
        }
        m_weight += w;
        return true;
    }

    ref<VectorBlock> clone() const {
        ref<VectorBlock> clone = new VectorBlock(m_dimensions);
        copyTo(clone);
        return clone;
    }

    void copyTo(VectorBlock *copy) const {
        copy->m_dimensions = m_dimensions;
        memcpy(&copy->m_data[0], &m_data[0], m_dimensions * sizeof(Float));
        copy->m_weight = m_weight;
    }
     
    // ======================================================================
	//! @{ \name Implementation of the WorkResult interface
	// ======================================================================

	void load(Stream *stream);
	void save(Stream *stream) const;
	std::string toString() const;

	//! @}
	// ======================================================================

    MTS_DECLARE_CLASS()
protected:
    virtual ~VectorBlock() {}
protected:
    Point2i m_offset;
    int m_dimensions;
    std::vector<Float> m_data;
    double m_weight;
};


MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_VECTORBLOCK_H_ */
