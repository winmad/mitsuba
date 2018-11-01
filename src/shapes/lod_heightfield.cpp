#include <mitsuba/render/shape.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sensor.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/timer.h>

MTS_NAMESPACE_BEGIN

class LodHeightField : public Shape {
public:
	LodHeightField(const Properties &props) : Shape(props) {
		m_levels = props.getInteger("levels", m_levels);
		m_effReso = Vector2i(
			props.getInteger("effResoX", 0),
			props.getInteger("effResoY", 0));
		
		m_footprintSizes.resize(m_levels);
		for (int i = 0; i < m_levels; i++) {
			std::ostringstream oss;
			oss << "footprintSize_" << i;
			m_footprintSizes[i] = props.getInteger(oss.str(), 0);
		}

		m_lods.resize(m_levels);
	}

	LodHeightField(Stream *stream, InstanceManager *manager)
		: Shape(stream, manager) {
		Log(EWarn, "Constructor from stream not implemented!");
		NotImplementedError("LodHeightField");
	}

	~LodHeightField() {
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Log(EWarn, "Serialize to stream not implemented!");
		NotImplementedError("serialize");
	}

	AABB getAABB() const {
		AABB result = m_lods[0]->getAABB();
		for (int i = 1; i < m_levels; i++) {
			result.expandBy(m_lods[i]->getAABB());
		}
		return result;
	}

	Float getSurfaceArea() const {
		return m_lods[m_levels - 1]->getSurfaceArea();
	}

	size_t getPrimitiveCount() const {
		return 1;
	}

	size_t getEffectivePrimitiveCount() const {
		return m_lods[m_levels - 1]->getEffectivePrimitiveCount();
	}

	inline Float getContinuousLevel(const Ray &ray) const {
		Float texels = ray.uvFootprint * m_effTexReso;
		if (texels >= m_footprintSizes[m_levels - 1])
			return m_levels - 1;
		if (texels < m_footprintSizes[0])
			return 0.0;
		for (int i = m_levels - 2; i >= 0; i--) {
			if (texels >= m_footprintSizes[i]) {
				Float x = (Float)(texels - m_footprintSizes[i]) / (m_footprintSizes[i + 1] - m_footprintSizes[i]);
				return i + x;
			}
		}
	}

	inline int chooseLevel(const Ray &ray, Float *lodLevel=NULL) const {
		Float cLevel = getContinuousLevel(ray);
		int res = math::clamp(math::floorToInt(cLevel) + ray.levelOffset, 0, m_levels - 1);
		if (lodLevel)
			*lodLevel = cLevel;
		return res;
	}

	void fillIntersectionRecord(const Ray &ray,
		const void *tmp, Intersection &its) const {
		Float lodLevel = 0.0;
		int levelIdx = chooseLevel(ray, &lodLevel);
		m_lods[levelIdx]->fillIntersectionRecord(ray, tmp, its);
	
		its.effTexReso = m_effTexReso;	
		its.lodLevel = lodLevel;
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt, Float &t, void *temp) const {
		int levelIdx = chooseLevel(ray);
		return m_lods[levelIdx]->rayIntersect(ray, mint, maxt, t, temp);
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
		Float t;
		return rayIntersect(ray, mint, maxt, t, NULL);
	}

	void getNormalDerivative(const Intersection &its, 
		Vector &dndu, Vector &dndv, bool shadingFrame) const {
		Log(EWarn, "getNormalDerivative not implemented!");
		NotImplementedError("getNormalDerivative");
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		for (int i = 0; i < m_levels; i++) {
			std::ostringstream oss;
			oss << "lod_" << i;
			if (name == oss.str()) {
				Log(EInfo, "%s", name.c_str());
				m_lods[i] = static_cast<Shape *>(child);
				return;
			}
		}
		Shape::addChild(name, child);
	}

	Shape *getElement(int i) {
		if (i < 0 || i >= m_levels)
			return NULL;
		return m_lods[i].get();
	}

	void configure() {
		Shape::configure();
		Log(EInfo, "%s", toString().c_str());

		m_effTexReso = std::min(m_effReso.x, m_effReso.y);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "LodHeightField[" << endl
			<< " levels = " << m_levels << "," << endl
			<< " effReso = (" << m_effReso.x << ", " << m_effReso.y << ")" << endl
			<< " footprintSizes:";
		
		for (int i = 0; i < m_levels; i++) {
			oss << " " << m_footprintSizes[i];
		}
		oss << endl << "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	int m_levels;
	Vector2i m_effReso;
	int m_effTexReso;
	std::vector<int> m_footprintSizes;
	ref_vector<Shape> m_lods;
};

MTS_IMPLEMENT_CLASS_S(LodHeightField, false, Shape)
MTS_EXPORT_PLUGIN(LodHeightField, "LoD height field intersection shape");
MTS_NAMESPACE_END