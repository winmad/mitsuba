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

#define MTS_QTREE_MAXDEPTH  50
#define MTS_QTREE_FASTSTART 1

//#define TILED_HEIGHTFIELD_DEBUG

MTS_NAMESPACE_BEGIN

static StatsCounter numTraversals("Height field", "Traversal operations per query", EAverage);

namespace {
	/// Find the smallest t >= 0 such that a*t + b is a multiple of c
	inline Float nextMultiple(Float a, Float b, Float c) {
		Float tmp = b / c,
			rounded = (a > 0 ? std::ceil(tmp) : std::floor(tmp)) * c,
			diff = rounded - b;

		if (diff == 0)
			diff = math::signum(a) * c;

		return diff / a;
	}

	/// Temporary storage for patch-ray intersections
	struct PatchIntersectionRecord {
		Point p;
		int x, y;
		int blockX, blockY;
	};

	/// Stack entry for recursive quadtree traversal
	struct StackEntry {
		int level, x, y;
	};
};

class TiledHeightfield : public Shape {
public:
	TiledHeightfield(const Properties &props) : Shape(props), m_data(NULL), m_normals(NULL), m_minmax(NULL) {
		m_sizeHint = Vector2i(
			props.getInteger("width", -1),
			props.getInteger("height", -1)
			);

		m_objectToWorld = props.getTransform("toWorld", Transform());
		m_shadingNormals = props.getBoolean("shadingNormals", false);
		m_flipNormals = props.getBoolean("flipNormals", false);
		m_scale = props.getFloat("scale", 1);
		m_tileX = props.getFloat("tileX", 1);
		m_tileY = props.getFloat("tileY", 1);

		m_filename = props.getString("filename", "");
		if (!m_filename.empty())
			m_filename = Thread::getThread()->getFileResolver()->resolve(m_filename);
	}

	TiledHeightfield(Stream *stream, InstanceManager *manager)
		: Shape(stream, manager), m_data(NULL), m_normals(NULL), m_minmax(NULL) {

		m_objectToWorld = Transform(stream);
		m_shadingNormals = stream->readBool();
		m_flipNormals = stream->readBool();
		m_scale = stream->readFloat();
		m_tileX = stream->readFloat();
		m_tileY = stream->readFloat();
		m_filename = stream->readString();
		m_dataSize = Vector2i(stream);
		size_t size = (size_t)m_dataSize.x * (size_t)m_dataSize.y;
		m_data = (Float *)allocAligned(size * sizeof(Float));
		stream->readFloatArray(m_data, size);
		configure();
	}

	~TiledHeightfield() {
		if (m_data)
			freeAligned(m_data);
		if (m_minmax) {
			for (int i = 0; i<m_levelCount; ++i)
				freeAligned(m_minmax[i]);
			delete[] m_minmax;
			delete[] m_levelSize;
			delete[] m_numChildren;
			delete[] m_blockSize;
			delete[] m_blockSizeF;
		}
		if (m_normals)
			freeAligned(m_normals);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
		m_objectToWorld.serialize(stream);
		stream->writeBool(m_shadingNormals);
		stream->writeBool(m_flipNormals);
		stream->writeFloat(m_scale);
		stream->writeFloat(m_tileX);
		stream->writeFloat(m_tileY);
		stream->writeString(m_filename.string());
		m_dataSize.serialize(stream);
		stream->writeFloatArray(m_data, (size_t)m_dataSize.x * (size_t)m_dataSize.y);
	}

	AABB getAABB() const {
		AABB result;
		for (int i = 0; i<8; ++i)
			result.expandBy(m_objectToWorld(m_totAABB.getCorner(i)));

		return result;
	}

	Float getSurfaceArea() const {
		return m_surfaceArea; /// XXX transformed surface area? ...
	}

	size_t getPrimitiveCount() const {
		return 1;
	}

	size_t getEffectivePrimitiveCount() const {
		return (size_t)m_levelSize[0].x * (size_t)m_levelSize[0].y;
	}

	inline static int signumToInt(Float value) {
		if (value < 0)
			return -1;
		else if (value > 0)
			return 1;
		else
			return 0;
	}

	bool rayIntersectLocal(const Ray &_ray, Float mint, Float maxt, Float &t, 
		Float blockX, Float blockY, void *tmp) const {
		StackEntry stack[MTS_QTREE_MAXDEPTH];

		Ray ray(_ray);

		/* Ray length to cross a single cell along the X or Y axis */
		Float tDeltaXSingle = std::abs(ray.dRcp.x),
			tDeltaYSingle = std::abs(ray.dRcp.y);

		/* Cell coordinate increments for steps along the ray */
		int iDeltaX = signumToInt(ray.d.x),
			iDeltaY = signumToInt(ray.d.y);

		int stackIdx = 0;

#if MTS_QTREE_FASTSTART
		/* If the entire ray is restricted to a subtree of the quadtree,
		directly start the traversal from the there instead of the root
		node. This can save some unnecessary work. */
		{
			Point enterPt, exitPt;
			Float nearT = mint, farT = maxt;
			if (!m_dataAABB.rayIntersect(ray, nearT, farT, enterPt, exitPt))
				return false;

			/* Determine minima and maxima in integer coordinates (round down!) */
			int minX = (int)std::min(enterPt.x, exitPt.x),
				maxX = (int)std::max(enterPt.x, exitPt.x),
				minY = (int)std::min(enterPt.y, exitPt.y),
				maxY = (int)std::max(enterPt.y, exitPt.y);

			/* Determine quadtree level */
			int level = math::clamp(1 + math::log2i(
				std::max((uint32_t)(minX ^ maxX), (uint32_t)(minY ^ maxY))),
				0, m_levelCount - 1);

			/* Compute X and Y coordinates at that level */
			const Vector2i &blockSize = m_blockSize[level];
			int x = math::clamp(minX / blockSize.x, 0, m_levelSize[level].x - 1),
				y = math::clamp(minY / blockSize.y, 0, m_levelSize[level].y - 1);

			stack[stackIdx].level = level;
			stack[stackIdx].x = x;
			stack[stackIdx].y = y;
		}
#else
		/* Start traversal from the root node of the quadtree */
		stack[stackIdx].level = m_levelCount - 1;
		stack[stackIdx].x = 0;
		stack[stackIdx].y = 0;
#endif

		numTraversals.incrementBase();

		size_t nTraversals = 0;
		while (stackIdx >= 0) {
			++nTraversals;

			/* Pop a node from the stack and compute its bounding box */
			StackEntry entry = stack[stackIdx--];
			const Interval &interval = m_minmax[entry.level][
				entry.x + entry.y * m_levelSize[entry.level].x];
			const Vector2 &blockSize = m_blockSizeF[entry.level];
			AABB aabb(
				Point3(0, 0, interval.min),
				Point3(blockSize.x, blockSize.y, interval.max)
				);

			/* Intersect the ray against the bounding box, in local coordinates */
			Ray localRay(Point(ray.o.x - entry.x*blockSize.x,
				ray.o.y - entry.y*blockSize.y, ray.o.z), ray.d, 0);
			Float nearT = mint, farT = maxt;
			Point enterPt, exitPt;

			if (!aabb.rayIntersect(localRay, nearT, farT, enterPt, exitPt)) {
				/* The bounding box was not intersected -- skip */
				continue;
			}

			Float tMax = farT - nearT;

			if (entry.level > 0) {
				/* Inner node -- push child nodes in 2D DDA order */
				const Vector2i &numChildren = m_numChildren[entry.level];
				const Vector2 &subBlockSize = m_blockSizeF[--entry.level];
				entry.x *= numChildren.x; entry.y *= numChildren.y;

				int x = (exitPt.x >= subBlockSize.x) ? numChildren.x - 1 : 0;
				int y = (exitPt.y >= subBlockSize.y) ? numChildren.y - 1 : 0;

				Float tDeltaX = tDeltaXSingle * subBlockSize.x,
					tDeltaY = tDeltaYSingle * subBlockSize.y,
					tNextX = nextMultiple(-ray.d.x, exitPt.x, subBlockSize.x),
					tNextY = nextMultiple(-ray.d.y, exitPt.y, subBlockSize.y),
					t = 0;

				while ((uint32_t)x < (uint32_t)numChildren.x &&
					(uint32_t)y < (uint32_t)numChildren.y && t <= tMax) {
					stack[++stackIdx].level = entry.level;
					stack[stackIdx].x = entry.x + x;
					stack[stackIdx].y = entry.y + y;

					if (tNextX < tNextY) {
						t = tNextX;
						tNextX += tDeltaX;
						x -= iDeltaX;
					}
					else {
						t = tNextY;
						tNextY += tDeltaY;
						y -= iDeltaY;
					}
				}
			}
			else {
				/* Intersect the ray against a bilinear patch */
				Float
					f00 = m_data[entry.y * m_dataSize.x + entry.x],
					f01 = m_data[((entry.y + 1) % m_dataSize.y) * m_dataSize.x + entry.x],
					f10 = m_data[entry.y * m_dataSize.x + (entry.x + 1) % m_dataSize.x],
					f11 = m_data[((entry.y + 1) % m_dataSize.y) * m_dataSize.x + (entry.x + 1) % m_dataSize.x];

				Float A = ray.d.x * ray.d.y * (f00 - f01 - f10 + f11);
				Float B = ray.d.y * (f01 - f00 + enterPt.x * (f00 - f01 - f10 + f11))
					+ ray.d.x * (f10 - f00 + enterPt.y * (f00 - f01 - f10 + f11))
					- ray.d.z;
				Float C = (enterPt.x - 1) * (enterPt.y - 1) * f00
					+ enterPt.y * f01 + enterPt.x * (f10 - enterPt.y * (f01 + f10 - f11))
					- enterPt.z;

				Float t0, t1;
				if (!solveQuadratic(A, B, C, t0, t1))
					continue;

				Float min = std::max(-Epsilon, mint - nearT);
				Float max = std::min(tMax + Epsilon, maxt - nearT);

				if (t0 >= min && t0 <= max)
					t = t0;
				else if (t1 >= min && t1 <= max)
					t = t1;
				else
					continue;

				if (tmp) {
					PatchIntersectionRecord &temp = *((PatchIntersectionRecord *)tmp);
					Point pLocal = enterPt + ray.d * t;
					temp.x = entry.x;
					temp.y = entry.y;
					temp.p = pLocal;
					temp.blockX = blockX;
					temp.blockY = blockY;
					t += nearT;
				}
				numTraversals += nTraversals;
				return true;
			}
		}

		numTraversals += nTraversals;
		return false;
	}

	void fillIntersectionRecord(const Ray &ray,
		const void *tmp, Intersection &its) const {
		PatchIntersectionRecord &temp = *((PatchIntersectionRecord *)tmp);

		int x = temp.x, y = temp.y, width = m_dataSize.x;
		Float
			f00 = m_data[y     * width + x],
			f01 = m_data[((y + 1) % m_dataSize.y) * width + x],
			f10 = m_data[y     * width + (x + 1) % m_dataSize.x],
			f11 = m_data[((y + 1) % m_dataSize.y) * width + (x + 1) % m_dataSize.x];

		Point pLocal(temp.p.x + temp.x + temp.blockX * m_dataSize.x, 
			temp.p.y + temp.y + temp.blockY * m_dataSize.y, temp.p.z);
		
		its.uv = Point2(pLocal.x * m_invTotSize.x, pLocal.y * m_invTotSize.y);
		its.p = m_objectToWorld(pLocal);
		its.dpdu = m_objectToWorld(Vector(1, 0,
			(1.0f - temp.p.y) * (f10 - f00) + temp.p.y * (f11 - f01)) * (m_levelSize0f.x * m_tileX));
		its.dpdv = m_objectToWorld(Vector(0, 1,
			(1.0f - temp.p.x) * (f01 - f00) + temp.p.x * (f11 - f10)) * (m_levelSize0f.y * m_tileY));

		its.geoFrame.s = normalize(its.dpdu);
		its.geoFrame.t = normalize(its.dpdv - dot(its.dpdv, its.geoFrame.s) * its.geoFrame.s);
		its.geoFrame.n = cross(its.geoFrame.s, its.geoFrame.t);

		if (m_shadingNormals) {
			const Normal
				&n00 = m_normals[y     * (width + 1) + x],
				&n01 = m_normals[(y + 1) * (width + 1) + x],
				&n10 = m_normals[y     * (width + 1) + (x + 1)],
				&n11 = m_normals[(y + 1) * (width + 1) + (x + 1)];

			its.shFrame.n = normalize(m_objectToWorld(Normal(
				(1 - temp.p.x) * ((1 - temp.p.y) * n00 + temp.p.y * n01)
				+ temp.p.x  * ((1 - temp.p.y) * n10 + temp.p.y * n11))));
		}
		else {
			its.shFrame.n = its.geoFrame.n;
		}

		if (m_flipNormals) {
			its.shFrame.n *= -1;
			its.geoFrame.n *= -1;
		}

		its.shape = this;
		its.hasUVPartials = false;
		its.instance = NULL;
		its.time = ray.time;
		its.primIndex = x + y * width;
		its.blockIndex = temp.blockX + temp.blockY * m_tileX;

		its.baseFrame = m_baseFrame;
	}

	bool rayIntersect(const Ray &_ray, Float mint, Float maxt, Float &t, void *temp) const {
		Ray ray;
		Point posWorld = _ray.o;
		m_worldToObject(_ray, ray);

		// check if the ray origin is below the heightfield
		/*
		if (m_basePlaneAABB.contains(Point2(ray.o.x, ray.o.y))) {
			Float h = getHeight(_ray.o);
			if (ray.o.z > -Epsilon && ray.o.z < h + Epsilon)
				return false;
		}
		*/

#ifdef TILED_HEIGHTFIELD_DEBUG
 		Log(EInfo, "------------------------------------");
 		Log(EInfo, "ray world: o = (%.3f, %.3f, %.3f), d = (%.3f, %.3f, %.3f)",
 			_ray.o.x, _ray.o.y, _ray.o.z, _ray.d.x, _ray.d.y, _ray.d.z);
 		Log(EInfo, "ray local: o = (%.3f, %.3f, %.3f), d = (%.3f, %.3f, %.3f)",
 			ray.o.x, ray.o.y, ray.o.z, ray.d.x, ray.d.y, ray.d.z);
		Log(EInfo, "mint = %.6f, maxt = %.6f", mint, maxt);
#endif

		Point enterPt, exitPt;
		Float nearT = mint, farT = maxt;
		if (!m_totAABB.rayIntersect(ray, nearT, farT, enterPt, exitPt))
			return false;

		t = 0;
		Float eps = Epsilon * 1.0f;
		if (nearT > Epsilon) {
			posWorld += _ray.d * (nearT + eps);
			mint -= nearT + eps;
			maxt -= nearT + eps;
			t += nearT + eps;
		}

		while (true) {
			Point posLocal = m_worldToObject.transformAffine(posWorld);
			Float tLocal;
			if (!m_totAABB.contains(posLocal))
				break;

			int blockX = math::floorToInt(posLocal.x / m_dataSize.x);
			int blockY = math::floorToInt(posLocal.y / m_dataSize.y);
			ray.o.x = posLocal.x - m_dataSize.x * blockX;
			ray.o.y = posLocal.y - m_dataSize.y * blockY;
			ray.o.z = posLocal.z;

#ifdef TILED_HEIGHTFIELD_DEBUG
 			Log(EInfo, "===== local block intersect =====");
 			Log(EInfo, "block: (%d, %d)", blockX, blockY);
 			Log(EInfo, "ray in local block: o = (%.3f, %.3f, %.3f), d = (%.3f, %.3f, %.3f)",
 				ray.o.x, ray.o.y, ray.o.z, ray.d.x, ray.d.y, ray.d.z);
			Log(EInfo, "mint = %.3f, maxt = %.3f", mint, maxt);
#endif 

			bool flag = rayIntersectLocal(ray, mint, maxt, tLocal, blockX, blockY, temp);

			if (flag) {
#ifdef TILED_HEIGHTFIELD_DEBUG
 				if (temp) {
 					PatchIntersectionRecord &tmp = *((PatchIntersectionRecord *)temp);
 					Log(EInfo, "texel = (%d, %d), uv = (%.3f, %.3f)",
 						tmp.x, tmp.y, tmp.p.x, tmp.p.y);
				}
#endif
				t += tLocal;
				return true;
			}

			AABB xyAABB(m_dataAABB);
			nearT = Epsilon; farT = std::numeric_limits<Float>::infinity();
			xyAABB.rayIntersect(ray, nearT, farT);
#ifdef TILED_HEIGHTFIELD_DEBUG
			Log(EInfo, "%.6f, %.6f", nearT, farT);
#endif
			Assert(nearT < Epsilon);
			posWorld += _ray.d * (farT + eps);
			mint -= farT + eps;
			maxt -= farT + eps;
			t += farT + eps;
		}
		
		t = std::numeric_limits<Float>::infinity();
		return false;
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
		Float t;
		return rayIntersect(ray, mint, maxt, t, NULL);
	}

	void getNormalDerivative(const Intersection &its,
		Vector &dndu, Vector &dndv, bool shadingFrame) const {
		int width = m_dataSize.x,
			blockX = its.blockIndex % m_tileX,
			blockY = its.blockIndex / m_tileX,
			x = its.primIndex % width,
			y = its.primIndex / width;

		Float u = its.uv.x * (m_levelSize0f.x * m_tileX) - (blockX * m_levelSize0f.x + x);
		Float v = its.uv.y * (m_levelSize0f.y * m_tileY) - (blockY * m_levelSize0f.y + y);

		Normal normal;
		if (shadingFrame && m_shadingNormals) {
			/* Derivatives for bilinear patch with interpolated shading normals */
			const Normal
				&n00 = m_normals[y     * (width + 1) + x],
				&n01 = m_normals[(y + 1) * (width + 1) + x],
				&n10 = m_normals[y     * (width + 1) + (x + 1)],
				&n11 = m_normals[(y + 1) * (width + 1) + (x + 1)];

			normal = m_objectToWorld(Normal(
				(1 - u) * ((1 - v) * n00 + v * n01)
				+ u  * ((1 - v) * n10 + v * n11)));

			dndu = m_objectToWorld(Normal((1.0f - v) * (n10 - n00) + v * (n11 - n01))) * (m_levelSize0f.x * m_tileX);
			dndv = m_objectToWorld(Normal((1.0f - u) * (n01 - n00) + u * (n11 - n10))) * (m_levelSize0f.y * m_tileY);
		}
		else {
			/* Derivatives for bilinear patch with geometric normals */
			Float
				f00 = m_data[y     * width + x],
				f01 = m_data[((y + 1) % m_dataSize.y) * width + x],
				f10 = m_data[y     * width + (x + 1) % m_dataSize.x],
				f11 = m_data[((y + 1) % m_dataSize.y) * width + (x + 1) % m_dataSize.x];

			normal = m_objectToWorld(
				Normal(f00 - f10 + (f01 + f10 - f00 - f11)*v,
				f00 - f01 + (f01 + f10 - f00 - f11)*u, 1));

			dndu = m_objectToWorld(Normal(0, f01 + f10 - f00 - f11, 0)) * (m_levelSize0f.x * m_tileX);
			dndv = m_objectToWorld(Normal(f01 + f10 - f00 - f11, 0, 0)) * (m_levelSize0f.y * m_tileY);
		}

		/* Account for normalization */
		Float invLength = 1 / normal.length();

		normal *= invLength;
		dndu *= invLength;
		dndv *= invLength;

		dndu -= dot(normal, dndu) * normal;
		dndv -= dot(normal, dndv) * normal;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();
		if (cClass->derivesFrom(Texture::m_theClass)) {
			if (m_data != NULL)
				Log(EError, "Attempted to attach multiple textures to a height field shape!");

			m_bitmap = static_cast<Texture *>(child)->getBitmap(m_sizeHint);
		}
		else if (cClass->derivesFrom(ReconstructionFilter::m_theClass)) {
			if (m_rfilter != NULL)
				Log(EError, "Attempted to attach multiple reconstruction filters to a height field shape!");

			m_rfilter = static_cast<ReconstructionFilter *>(child);
		}
		else {
			Shape::addChild(name, child);
		}
	}

	void configure() {
		Shape::configure();

		if (m_minmax)
			return;

		if (!m_filename.empty()) {
			if (m_bitmap.get())
				Log(EError, "Cannot specify a file name and a nested texture at the same time!");
			ref<FileStream> fs = new FileStream(m_filename, FileStream::EReadOnly);
			m_bitmap = new Bitmap(Bitmap::EAuto, fs);
		}
		else if (!m_bitmap.get()) {
			Log(EError, "A height field texture must be specified (either as a nested texture, or using the 'filename' parameter)");
		}

		m_dataSize = m_bitmap->getSize();
		if (m_dataSize.x < 2) m_dataSize.x = 2;
		if (m_dataSize.y < 2) m_dataSize.y = 2;
		if (!math::isPowerOfTwo(m_dataSize.x)) m_dataSize.x = (int)math::roundToPowerOfTwo((uint32_t)m_dataSize.x);
		if (!math::isPowerOfTwo(m_dataSize.y)) m_dataSize.y = (int)math::roundToPowerOfTwo((uint32_t)m_dataSize.y);

		if (m_bitmap->getSize() != m_dataSize) {
			m_bitmap = m_bitmap->convert(Bitmap::ELuminance, Bitmap::EFloat);

			Log(EInfo, "Resampling heightfield texture from %ix%i to %ix%i ..",
				m_bitmap->getWidth(), m_bitmap->getHeight(), m_dataSize.x, m_dataSize.y);

			m_bitmap = m_bitmap->resample(m_rfilter, ReconstructionFilter::EClamp,
				ReconstructionFilter::EClamp, m_dataSize,
				-std::numeric_limits<Float>::infinity(),
				std::numeric_limits<Float>::infinity());
		}

		m_totDataSize = Vector2i(m_dataSize.x * m_tileX, m_dataSize.y * m_tileY);

		size_t size = (size_t)m_dataSize.x * (size_t)m_dataSize.y * sizeof(Float);
		m_data = (Float *)allocAligned(size);
		m_bitmap->convert(m_data, Bitmap::ELuminance, Bitmap::EFloat, 1.0f, m_scale);

		m_objectToWorld = m_objectToWorld * Transform::translate(Vector(-1, -1, 0)) * 
			Transform::scale(Vector((Float)1 / m_dataSize.x, (Float)1 / m_dataSize.y, 1)) *
			Transform::scale(Vector((Float)2 / m_tileX, (Float)2 / m_tileY, 1));
		m_worldToObject = m_objectToWorld.inverse();

		m_bitmap = NULL;

		size_t storageSize = (size_t)m_dataSize.x * (size_t)m_dataSize.y * sizeof(Float);
		Log(EInfo, "Building acceleration data structure for %ix%i height field ..", m_dataSize.x, m_dataSize.y);

		ref<Timer> timer = new Timer();
		m_levelCount = (int)std::max(math::log2i((uint32_t)m_dataSize.x), math::log2i((uint32_t)m_dataSize.y)) + 1;

		m_levelSize = new Vector2i[m_levelCount];
		m_numChildren = new Vector2i[m_levelCount];
		m_blockSize = new Vector2i[m_levelCount];
		m_blockSizeF = new Vector2[m_levelCount];
		m_minmax = new Interval*[m_levelCount];

		m_levelSize[0] = Vector2i(m_dataSize.x, m_dataSize.y);
		m_levelSize0f = Vector2(m_levelSize[0]);
		m_blockSize[0] = Vector2i(1, 1);
		m_blockSizeF[0] = Vector2(1, 1);
		m_invSize = Vector2((Float)1 / m_levelSize[0].x, (Float)1 / m_levelSize[0].y);
		m_invTotSize = Vector2(m_invSize.x / m_tileX, m_invSize.y / m_tileY);
		m_surfaceArea = 0;

		size = (size_t)m_levelSize[0].x * (size_t)m_levelSize[0].y * sizeof(Interval);
		m_minmax[0] = (Interval *)allocAligned(size);
		storageSize += size;

		/* Build the lowest MIP layer directly from the heightfield data */
		Interval *bounds = m_minmax[0];
		for (int y = 0; y<m_levelSize[0].y; ++y) {
			for (int x = 0; x<m_levelSize[0].x; ++x) {
				Float f00 = m_data[y * m_dataSize.x + x];
				Float f10 = m_data[y * m_dataSize.x + (x + 1) % m_dataSize.x];
				Float f01 = m_data[((y + 1) % m_dataSize.y) * m_dataSize.x + x];
				Float f11 = m_data[((y + 1) % m_dataSize.y) * m_dataSize.x + (x + 1) % m_dataSize.x];
				Float fmin = std::min(std::min(f00, f01), std::min(f10, f11));
				Float fmax = std::max(std::max(f00, f01), std::max(f10, f11));
				*bounds++ = Interval(fmin, fmax);

				/* Estimate the total surface area (this is approximate) */
				Float diff0 = f01 - f10, diff1 = f00 - f11;
				m_surfaceArea += std::sqrt(1.0f + .5f * (diff0*diff0 + diff1*diff1));
			}
		}

		/* Propagate height bounds upwards to the other layers */
		for (int level = 1; level<m_levelCount; ++level) {
			Vector2i &cur = m_levelSize[level],
				&prev = m_levelSize[level - 1];

			/* Calculate size of this layer */
			cur.x = prev.x > 1 ? (prev.x / 2) : 1;
			cur.y = prev.y > 1 ? (prev.y / 2) : 1;

			m_numChildren[level].x = prev.x > 1 ? 2 : 1;
			m_numChildren[level].y = prev.y > 1 ? 2 : 1;
			m_blockSize[level] = Vector2i(
				m_levelSize[0].x / cur.x,
				m_levelSize[0].y / cur.y
				);
			m_blockSizeF[level] = Vector2(m_blockSize[level]);

			/* Allocate memory for interval data */
			Interval *prevBounds = m_minmax[level - 1], *curBounds;
			size_t size = (size_t)cur.x * (size_t)cur.y * sizeof(Interval);
			m_minmax[level] = curBounds = (Interval *)allocAligned(size);
			storageSize += size;

			/* Build by querying the previous layer */
			for (int y = 0; y<cur.y; ++y) {
				int y0 = std::min(2 * y, prev.y - 1),
					y1 = std::min(2 * y + 1, prev.y - 1);
				for (int x = 0; x<cur.x; ++x) {
					int x0 = std::min(2 * x, prev.x - 1),
						x1 = std::min(2 * x + 1, prev.x - 1);
					const Interval &f00 = prevBounds[y0 * prev.x + x0], &f01 = prevBounds[y0 * prev.x + x1];
					const Interval &f10 = prevBounds[y1 * prev.x + x0], &f11 = prevBounds[y1 * prev.x + x1];
					Interval combined(f00);
					combined.expandBy(f01);
					combined.expandBy(f10);
					combined.expandBy(f11);
					*curBounds++ = combined;
				}
			}
		}

		if (m_shadingNormals) {
			Log(EInfo, "Precomputing shading normals ..");
			size_t size = (size_t)(m_dataSize.x + 1) * (size_t)(m_dataSize.y + 1) * sizeof(Normal);
			m_normals = (Normal *)allocAligned(size);
			memset(m_normals, 0, size);
			storageSize += size;

			for (int offset = 0; offset<2; ++offset) {
#if defined(MTS_OPENMP)
#pragma omp parallel for
#endif
				for (int y = offset; y<m_levelSize[0].y; y += 2) {
					for (int x = 0; x<m_levelSize[0].x; ++x) {
						Float f00 = m_data[y * m_dataSize.x + x];
						Float f10 = m_data[y * m_dataSize.x + (x + 1) % m_dataSize.x];
						Float f01 = m_data[((y + 1) % m_dataSize.y) * m_dataSize.x + x];
						Float f11 = m_data[((y + 1) % m_dataSize.y) * m_dataSize.x + (x + 1) % m_dataSize.x];

						m_normals[y       * (m_dataSize.x + 1) + x] += normalize(Normal(f00 - f10, f00 - f01, 1));
						m_normals[y       * (m_dataSize.x + 1) + (x + 1)] += normalize(Normal(f00 - f10, f10 - f11, 1));
						m_normals[(y + 1) * (m_dataSize.x + 1) + x] += normalize(Normal(f01 - f11, f00 - f01, 1));
						m_normals[(y + 1) * (m_dataSize.x + 1) + (x + 1)] += normalize(Normal(f01 - f11, f10 - f11, 1));
					}
				}
			}

#if defined(MTS_OPENMP)
#pragma omp parallel for
#endif
			for (int y = 0; y<=m_dataSize.y; ++y) {
				for (int x = 0; x<=m_dataSize.x; ++x) {
					Normal &normal = m_normals[x + y * (m_dataSize.x + 1)];
					normal /= normal.length();
				}
			}
		}

		Log(EInfo, "Done (took %i ms, uses %s of memory)", timer->getMilliseconds(),
			memString(storageSize).c_str());

		m_dataAABB = AABB(
			Point3(0, 0, m_minmax[m_levelCount - 1][0].min),
			Point3(m_levelSize0f.x, m_levelSize0f.y, m_minmax[m_levelCount - 1][0].max)
			);

		m_totAABB = AABB(
			Point3(0, 0, m_dataAABB.min.z),
			Point3(m_levelSize0f.x * m_tileX, m_levelSize0f.y * m_tileY, m_dataAABB.max.z));

		m_basePlaneAABB = AABB2(
			Point2(0, 0),
			Point2(m_levelSize0f.x * m_tileX, m_levelSize0f.y * m_tileY)
		);

		Vector baseX = normalize(m_objectToWorld(Vector(1.0, 0.0, 0.0)));
		Normal baseN = normalize(m_objectToWorld(Normal(0.0, 0.0, 1.0)));
		Vector baseY = cross(baseN, baseX);
		m_baseFrame = Frame(baseX, baseY, baseN);
	}

	ref<TriMesh> createTriMesh() {
		Vector2i size = m_dataSize;

		/* Limit the size of the mesh */
		while (size.x > 2048 && size.y > 2048) {
			size.x = std::max(size.x / 2, 2);
			size.y = std::max(size.y / 2, 2);
		}

		size_t numTris = 2 * (size_t)size.x * (size_t)size.y;
		size_t numVertices = (size_t)(size.x + 1) * (size_t)(size.y + 1);

		ref<TriMesh> mesh = new TriMesh("Height field approximation",
			numTris, numVertices, false, true, false, false, !m_shadingNormals);

		Point *vertices = mesh->getVertexPositions();
		Point2 *texcoords = mesh->getVertexTexcoords();
		Triangle *triangles = mesh->getTriangles();

		Float dx = (Float)1 / size.x;
		Float dy = (Float)1 / size.y;
		Float scaleX = (Float)m_dataSize.x / size.x;
		Float scaleY = (Float)m_dataSize.y / size.y;

		uint32_t vertexIdx = 0;
		for (int y = 0; y<=size.y; ++y) {
			int py = std::min((int)(scaleY * y), m_dataSize.y);
			for (int x = 0; x<=size.x; ++x) {
				int px = std::min((int)(scaleX * x), m_dataSize.x);
				texcoords[vertexIdx] = Point2(x*dx, y*dy);
				vertices[vertexIdx++] = m_objectToWorld(Point((Float)px, (Float)py,
					m_data[px % m_dataSize.x + (py % m_dataSize.y) * m_dataSize.x]));
			}
		}
		Assert(vertexIdx == numVertices);

		uint32_t width = size.x + 1;
		uint32_t triangleIdx = 0;
		for (int y = 1; y<=size.y; ++y) {
			for (int x = 0; x<=size.x - 1; ++x) {
				uint32_t nextx = x + 1;
				uint32_t idx0 = width*y + x;
				uint32_t idx1 = width*y + nextx;
				uint32_t idx2 = width*(y - 1) + x;
				uint32_t idx3 = width*(y - 1) + nextx;

				triangles[triangleIdx].idx[0] = idx0;
				triangles[triangleIdx].idx[1] = idx2;
				triangles[triangleIdx].idx[2] = idx1;
				triangleIdx++;
				triangles[triangleIdx].idx[0] = idx1;
				triangles[triangleIdx].idx[1] = idx2;
				triangles[triangleIdx].idx[2] = idx3;
				triangleIdx++;
			}
		}
		Assert(triangleIdx == numTris);
		mesh->copyAttachments(this);
		mesh->configure();

		return mesh.get();
	}

	Float getHeight(const Point &_p) const {
		Point p = m_worldToObject.transformAffine(_p);
		int x = math::floorToInt(p.x);
		int y = math::floorToInt(p.y);
		Float u = p.x - x;
		Float v = p.y - y;
		int blockX = math::floorToInt((Float)x / m_dataSize.x);
		int blockY = math::floorToInt((Float)y / m_dataSize.y);
		x -= m_dataSize.x * blockX;
		y -= m_dataSize.y * blockY;
		int width = m_dataSize.x;

		//Log(EInfo, "World (%.6f, %.6f, %.6f)", _p.x, _p.y, _p.z);
		//Log(EInfo, "Local (%.6f, %.6f, %.6f)", p.x, p.y, p.z);
		//Log(EInfo, "(%.6f, %.6f)", u, v);

		Float
			f00 = m_data[y     * width + x],
			f01 = m_data[((y + 1) % m_dataSize.y) * width + x],
			f10 = m_data[y     * width + (x + 1) % m_dataSize.x],
			f11 = m_data[((y + 1) % m_dataSize.y) * width + (x + 1) % m_dataSize.x];

		//Log(EInfo, "%.6f, %.6f, %.6f, %.6f", f00, f01, f10, f11);

		Float res = (1.0f - u) * (1.0f - v) * f00 + (1.0f - u) * v * f01 +
			u * (1.0f - v) * f10 + u * v * f11;
		return res;
	}

	Normal getNormal(const Point &_p) const {
		Point p = m_worldToObject.transformAffine(_p);
		int x = math::floorToInt(p.x);
		int y = math::floorToInt(p.y);
		Float u = p.x - x;
		Float v = p.y - y;
		int blockX = math::floorToInt((Float)x / m_dataSize.x);
		int blockY = math::floorToInt((Float)y / m_dataSize.y);
		x -= m_dataSize.x * blockX;
		y -= m_dataSize.y * blockY;
		int width = m_dataSize.x;

		Float
			f00 = m_data[y     * width + x],
			f01 = m_data[((y + 1) % m_dataSize.y) * width + x],
			f10 = m_data[y     * width + (x + 1) % m_dataSize.x],
			f11 = m_data[((y + 1) % m_dataSize.y) * width + (x + 1) % m_dataSize.x];

		Normal normal = m_objectToWorld(
			Normal(f00 - f10 + (f01 + f10 - f00 - f11)*v,
			f00 - f01 + (f01 + f10 - f00 - f11)*u, 1));

// 		Log(EInfo, "World (%.6f, %.6f, %.6f)", _p.x, _p.y, _p.z);
// 		Log(EInfo, "Local (%.6f, %.6f, %.6f)", p.x, p.y, p.z);
// 		Log(EInfo, "(%.6f, %.6f)", u, v);
// 		Log(EInfo, "%.6f, %.6f, %.6f, %.6f", f00, f01, f10, f11);

		normal = normalize(normal);
		return normal;
	}

	void getPosAndNormal(const Point2 &uv, Point *pos, Normal *normal) const {
		Point p;
		p.x = uv.x * (m_levelSize0f.x * m_tileX);
		p.y = uv.y * (m_levelSize0f.y * m_tileY);
		int x = math::floorToInt(p.x);
		int y = math::floorToInt(p.y);
		Float u = p.x - x;
		Float v = p.y - y;
		int blockX = math::floorToInt((Float)x / m_dataSize.x);
		int blockY = math::floorToInt((Float)y / m_dataSize.y);
		x -= m_dataSize.x * blockX;
		y -= m_dataSize.y * blockY;
		int width = m_dataSize.x;

		Float
			f00 = m_data[y     * width + x],
			f01 = m_data[((y + 1) % m_dataSize.y) * width + x],
			f10 = m_data[y     * width + (x + 1) % m_dataSize.x],
			f11 = m_data[((y + 1) % m_dataSize.y) * width + (x + 1) % m_dataSize.x];

		p.z = (1.0f - u) * (1.0f - v) * f00 + (1.0f - u) * v * f01 +
			u * (1.0f - v) * f10 + u * v * f11;
		p = m_objectToWorld(p);

		Normal norm = m_objectToWorld(
			Normal(f00 - f10 + (f01 + f10 - f00 - f11)*v,
			f00 - f01 + (f01 + f10 - f00 - f11)*u, 1));
		norm = normalize(norm);

		if (pos != NULL) {
			*pos = p;
		}
		if (normal != NULL) {
			*normal = norm;
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "TiledHeightField[" << endl
			<< "  size = " << m_dataSize.toString() << "," << endl
			<< "  tileX = " << m_tileX << "," << endl
			<< "  tileY = " << m_tileY << "," << endl
			<< "  shadingNormals = " << m_shadingNormals << "," << endl
			<< "  flipNormals = " << m_flipNormals << "," << endl
			<< "  objectToWorld = " << indent(m_objectToWorld.toString()) << "," << endl
			<< "  aabb = " << indent(getAABB().toString()) << "," << endl
			<< "  bsdf = " << indent(m_bsdf.toString()) << "," << endl;
		if (isMediumTransition())
			oss << "  interiorMedium = " << indent(m_interiorMedium.toString()) << "," << endl
			<< "  exteriorMedium = " << indent(m_exteriorMedium.toString()) << "," << endl;
		oss << "  emitter = " << indent(m_emitter.toString()) << "," << endl
			<< "  sensor = " << indent(m_sensor.toString()) << "," << endl
			<< "  subsurface = " << indent(m_subsurface.toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<ReconstructionFilter> m_rfilter;
	ref<Bitmap> m_bitmap;
	Transform m_objectToWorld, m_worldToObject;
	Vector2i m_sizeHint;
	AABB m_dataAABB;
	bool m_shadingNormals;
	bool m_flipNormals;
	Float m_scale;
	fs::path m_filename;

	/* Height field data */
	Float *m_data;
	Normal *m_normals;
	Vector2i m_dataSize;
	Vector2 m_invSize;
	Float m_surfaceArea;

	/* Min-max quadtree data */
	int m_levelCount;
	Vector2i *m_levelSize;
	Vector2   m_levelSize0f;
	Vector2i *m_numChildren;
	Vector2i *m_blockSize;
	Vector2 *m_blockSizeF;
	Interval **m_minmax;

	int m_tileX, m_tileY;
	Vector2i m_totDataSize;
	Vector2 m_invTotSize;
	AABB m_totAABB;
	AABB2 m_basePlaneAABB;
	Frame m_baseFrame;
};

MTS_IMPLEMENT_CLASS_S(TiledHeightfield, false, Shape)
MTS_EXPORT_PLUGIN(TiledHeightfield, "Tiled height field intersection shape");
MTS_NAMESPACE_END
