#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/pmf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/mipmap.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

class TabulatedScaledBSDF : public BSDF {
public:
	TabulatedScaledBSDF(const Properties &props) : BSDF(props) {
		m_angularScaleFilename = props.getString("angularScaleFilename", "");
		m_spatialScaleFilename = props.getString("spatialScaleFilename", "");
		m_spatialInterp = props.getString("spatialInterp", "nearest");
		m_useMIS = props.getBoolean("useMIS", false);

		m_wiUseFullSphere = props.getBoolean("wiUseFullSphere", false);
		m_woUseFullSphere = props.getBoolean("woUseFullSphere", false);
		m_useGrad = props.getBoolean("useGrad", false);

		Float uvscale = props.getFloat("uvscale", 1.0f);
		m_uvScale = Vector2(
			props.getFloat("uscale", uvscale),
			props.getFloat("vscale", uvscale));
	}

	TabulatedScaledBSDF(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		m_angularScaleFilename = stream->readString();
		m_spatialScaleFilename = stream->readString();
		m_useMIS = stream->readBool();

		m_wiUseFullSphere = stream->readBool();
		m_woUseFullSphere = stream->readBool();
		m_useGrad = stream->readBool();

		m_uvScale = Vector2(stream);

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeString(m_angularScaleFilename);
		stream->writeString(m_spatialScaleFilename);
		stream->writeBool(m_useMIS);

		stream->writeBool(m_wiUseFullSphere);
		stream->writeBool(m_woUseFullSphere);
		stream->writeBool(m_useGrad);

		m_uvScale.serialize(stream);
	}

	void configure() {
		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide | ESpatiallyVarying);

		m_usesRayDifferentials = false;

		// load angular scales
		if (m_angularScaleFilename != "")
			m_angularScales = new Bitmap(fs::path(m_angularScaleFilename));
		else
			m_angularScales = NULL;

		if (m_angularScales != NULL) {
			m_lobeSize = m_angularScales->getSize();

			if (m_wiUseFullSphere)
				m_wiResolution = math::floorToInt(std::sqrt((Float)m_lobeSize.y * 0.5));
			else
				m_wiResolution = math::floorToInt(std::sqrt((Float)m_lobeSize.y));

			if (m_woUseFullSphere)
				m_woResolution = math::floorToInt(std::sqrt((Float)m_lobeSize.x * 0.5));
			else
				m_woResolution = math::floorToInt(std::sqrt((Float)m_lobeSize.x));

			Log(EInfo, "wiRes = %d, woRes = %d", m_wiResolution, m_woResolution);

			if (m_useMIS) {
				m_rowCondCdfs.resize(m_lobeSize.y, DiscreteDistribution(m_lobeSize.x));
				for (int r = 0; r < m_lobeSize.y; r++) {
					for (int c = 0; c < m_lobeSize.x; c++) {
						Spectrum spec = m_angularScales->getPixel(Point2i(c, r));
						m_rowCondCdfs[r].append(spec.average());
					}
					m_rowCondCdfs[r].normalize();
				}

				m_samplers.resize(233);
				m_samplers[0] = static_cast<Sampler *>(PluginManager::getInstance()->
					createObject(MTS_CLASS(Sampler), Properties("independent")));
				m_samplers[0]->configure();
				for (int i = 1; i < 233; i++) {
					m_samplers[i] = m_samplers[0]->clone();
				}
			}
		}

		// load spatial scales
		if (m_spatialScaleFilename != "")
			m_spatialScales = new Bitmap(fs::path(m_spatialScaleFilename));
		else
			m_spatialScales = NULL;

		if (m_spatialScales != NULL) {
			m_xyResolution = m_spatialScales->getSize();
			Log(EInfo, "xRes = %d, yRes = %d", m_xyResolution.x, m_xyResolution.y);
		}

		BSDF::configure();
	}

	Spectrum evalSpatialScale(const BSDFSamplingRecord &bRec) const {
		Point2 uv = transformUV(bRec.its.uv);
		
		if (m_spatialInterp == "nearest") {
			int xIdx = math::floorToInt(uv.x * m_xyResolution.x);
			int yIdx = math::floorToInt(uv.y * m_xyResolution.y);
			return m_spatialScales->getPixel(Point2i(xIdx, yIdx));
		} else if (m_spatialInterp == "bilinear") {
			Float x = uv.x * m_xyResolution.x;
			Float y = uv.y * m_xyResolution.y;

			// grid point at (x+0.5, y+0.5)
			int xIdx = math::floorToInt(x + 0.5);
			int yIdx = math::floorToInt(y + 0.5);

// 			Log(EInfo, "==================");
// 			Log(EInfo, "(x, y) = (%.6f, %.6f)", x, y);

			Spectrum res(0.f);

			for (int dy = -1; dy <= 0; dy++) {
				Float wv = 1.0 - std::abs(y - (yIdx + dy + 0.5));
				int yNow = (yIdx + m_xyResolution.y + dy) % m_xyResolution.y;
				for (int dx = -1; dx <= 0; dx++) {
					Float wu = 1.0 - std::abs(x - (xIdx + dx + 0.5));
					int xNow = (xIdx + m_xyResolution.x + dx) % m_xyResolution.x;
					Spectrum tmpValue = m_spatialScales->getPixel(Point2i(xNow, yNow));
					res += tmpValue * wu * wv;

					//Log(EInfo, "(xNow, yNow), wu, wv: (%d, %d), %.6f, %.6f", xNow, yNow, wu, wv);
				}
			}
			return res;
		}		
	}

	Spectrum evalScale(const BSDFSamplingRecord &bRec) const {
		Vector wiWorld = bRec.its.toWorld(bRec.wi);
		Vector wiMacro = bRec.its.baseFrame.toLocal(wiWorld);
		Vector woWorld = bRec.its.toWorld(bRec.wo);
		Vector woMacro = bRec.its.baseFrame.toLocal(woWorld);

		// Added Oct 23, 2018
// 		if (wiMacro.z <= 0 || woMacro.z <= 0)
// 			return Spectrum(0.0);

// 		if (woMacro.z <= 0)
// 			return Spectrum(0.0);

		int r1Offset = 0;
		if (wiMacro.z <= 0) {
			if (!m_wiUseFullSphere)
				return Spectrum(0.0);
			r1Offset = m_wiResolution;
			wiMacro.z = -wiMacro.z;
		}

		int r2Offset = 0;
		if (woMacro.z <= 0) {
			if (!m_woUseFullSphere)
				return Spectrum(0.0);
			r2Offset = m_woResolution;
			woMacro.z = -woMacro.z;
		}

		Point2 wiTex = warp::uniformHemisphereToSquareConcentric(wiMacro);
		Point2 woTex = warp::uniformHemisphereToSquareConcentric(woMacro);

		// piecewise bilinear
		int wiNumCells = m_wiResolution - 1;
		int woNumCells = m_woResolution - 1;

		int c1 = math::clamp(math::floorToInt(wiTex.x * wiNumCells), 0, wiNumCells - 1);
		int r1 = math::clamp(math::floorToInt(wiTex.y * wiNumCells), 0, wiNumCells - 1) + r1Offset;
		int c2 = math::clamp(math::floorToInt(woTex.x * woNumCells), 0, woNumCells - 1);
		int r2 = math::clamp(math::floorToInt(woTex.y * woNumCells), 0, woNumCells - 1) + r2Offset;

		Spectrum res(0.f);
		Float w(0.f);
		for (int dr1 = 0; dr1 < 2; dr1++) {
			Float v1 = wiTex.y * wiNumCells - (r1 % m_wiResolution);
			Float wv1 = std::abs(1.0 - dr1 - v1);
			
			for (int dc1 = 0; dc1 < 2; dc1++) {
				Float u1 = wiTex.x * wiNumCells - c1;
				Float wu1 = std::abs(1.0 - dc1 - u1);
				
				for (int dr2 = 0; dr2 < 2; dr2++) {
					Float v2 = woTex.y * woNumCells - (r2 % m_woResolution);
					Float wv2 = std::abs(1.0 - dr2 - v2);
					
					for (int dc2 = 0; dc2 < 2; dc2++) {
						Float u2 = woTex.x * woNumCells - c2;
						Float wu2 = std::abs(1.0 - dc2 - u2);

						int wiIdx = (r1 + dr1) * m_wiResolution + (c1 + dc1);
						int woIdx = (r2 + dr2) * m_woResolution + (c2 + dc2);

						Spectrum tmpValue = m_angularScales->getPixel(Point2i(woIdx, wiIdx));
						res += tmpValue * wv1 * wu1 * wv2 * wu2;
						w += wv1 * wu1 * wv2 * wu2;
					}
				}
			}
		}

		return res;
	}

	Spectrum evalScaleAndGrad(BSDFSamplingRecord &bRec) const {
		Vector wiWorld = bRec.its.toWorld(bRec.wi);
		Vector wiMacro = bRec.its.baseFrame.toLocal(wiWorld);
		Vector woWorld = bRec.its.toWorld(bRec.wo);
		Vector woMacro = bRec.its.baseFrame.toLocal(woWorld);

		int r1Offset = 0;
		if (wiMacro.z <= 0) {
			if (!m_wiUseFullSphere)
				return Spectrum(0.0);
			r1Offset = m_wiResolution;
			wiMacro.z = -wiMacro.z;
		}

		int r2Offset = 0;
		if (woMacro.z <= 0) {
			if (!m_woUseFullSphere)
				return Spectrum(0.0);
			r2Offset = m_woResolution;
			woMacro.z = -woMacro.z;
		}

		Point2 wiTex = warp::uniformHemisphereToSquareConcentric(wiMacro);
		Point2 woTex = warp::uniformHemisphereToSquareConcentric(woMacro);

		// piecewise bilinear
		int wiNumCells = m_wiResolution - 1;
		int woNumCells = m_woResolution - 1;

		int c1 = math::clamp(math::floorToInt(wiTex.x * wiNumCells), 0, wiNumCells - 1);
		int r1 = math::clamp(math::floorToInt(wiTex.y * wiNumCells), 0, wiNumCells - 1) + r1Offset;
		int c2 = math::clamp(math::floorToInt(woTex.x * woNumCells), 0, woNumCells - 1);
		int r2 = math::clamp(math::floorToInt(woTex.y * woNumCells), 0, woNumCells - 1) + r2Offset;

		Spectrum res(0.f);
		int k = 0;
		for (int dr1 = 0; dr1 < 2; dr1++) {
			Float v1 = wiTex.y * wiNumCells - (r1 % m_wiResolution);
			Float wv1 = std::abs(1.0 - dr1 - v1);
			
			for (int dc1 = 0; dc1 < 2; dc1++) {
				Float u1 = wiTex.x * wiNumCells - c1;
				Float wu1 = std::abs(1.0 - dc1 - u1);
				
				for (int dr2 = 0; dr2 < 2; dr2++) {
					Float v2 = woTex.y * woNumCells - (r2 % m_woResolution);
					Float wv2 = std::abs(1.0 - dr2 - v2);
					
					for (int dc2 = 0; dc2 < 2; dc2++) {
						Float u2 = woTex.x * woNumCells - c2;
						Float wu2 = std::abs(1.0 - dc2 - u2);

						int wiIdx = (r1 + dr1) * m_wiResolution + (c1 + dc1);
						int woIdx = (r2 + dr2) * m_woResolution + (c2 + dc2);

						Spectrum tmpValue = m_angularScales->getPixel(Point2i(woIdx, wiIdx));
						Float weight = wv1 * wu1 * wv2 * wu2;
						res += tmpValue * weight;

						bRec.smIndices[k] = wiIdx * (2 * m_woResolution * m_woResolution) + woIdx;
						bRec.smDiff[k] = Spectrum(weight);
						k++;
					}
				}
			}
		}

		for (int i = 0; i < 16; i++)
			for (int c = 0; c < 3; c++)
				bRec.smDiff[i][c] /= res[c];

		return res;
	}

	Spectrum evalScaleDisk(BSDFSamplingRecord &bRec) const {
		Vector wiWorld = bRec.its.toWorld(bRec.wi);
		Vector wiMacro = bRec.its.baseFrame.toLocal(wiWorld);
		Vector woWorld = bRec.its.toWorld(bRec.wo);
		Vector woMacro = bRec.its.baseFrame.toLocal(woWorld);

		// no interpolation
// 		int r1 = math::floorToInt((wiMacro.y + 1.0) * 0.5 * m_wiResolution);
// 		r1 = math::clamp(m_wiResolution - r1 - 1, 0, m_wiResolution - 1);
// 		int c1 = math::floorToInt((wiMacro.x + 1.0) * 0.5 * m_wiResolution);
// 		c1 = math::clamp(c1, 0, m_wiResolution - 1);
// 
// 		int r2 = math::floorToInt((woMacro.y + 1.0) * 0.5 * m_woResolution);
// 		r2 = math::clamp(m_woResolution - r2 - 1, 0, m_woResolution - 1);
// 		int c2 = math::floorToInt((woMacro.x + 1.0) * 0.5 * m_woResolution);
// 		c2 = math::clamp(c2, 0, m_woResolution - 1);
// 
// 		Spectrum result = m_angularScales->getPixel(Point2i(r2 * m_woResolution + c2, 
// 			r1 * m_wiResolution + c1));
// 		if (result[0] < 1e-4f) {
// 			result = Spectrum(1.0f);
// 		}
// 		//Log(EInfo, "%.6f, %.6f, %.6f", result[0], result[1], result[2]);
// 		return result;

		if (wiMacro.z <= 0 || woMacro.z <= 0)
			return Spectrum(0.0);

		int rWi = math::floorToInt((wiMacro.y + 1.0) * 0.5 * m_wiResolution + 0.5);
		rWi = math::clamp(m_wiResolution - rWi - 1, 0, m_wiResolution - 1);
		int cWi = math::floorToInt((wiMacro.x + 1.0) * 0.5 * m_wiResolution - 0.5);
		cWi = math::clamp(cWi, 0, m_wiResolution - 1);

		int rWo = math::floorToInt((woMacro.y + 1.0) * 0.5 * m_woResolution + 0.5);
		rWo = math::clamp(m_woResolution - rWo - 1, 0, m_woResolution - 1);
		int cWo = math::floorToInt((woMacro.x + 1.0) * 0.5 * m_woResolution - 0.5);
		cWo = math::clamp(cWo, 0, m_woResolution - 1);

		Spectrum res(0.f);
		Float weights = 0.f;

		// interpolation
		for (int drWi = 0; drWi < 2; drWi++) {
			int rWiNew = rWi + drWi;
// 			if (rWiNew < 0 || rWiNew >= m_wiResolution)
// 				continue;
			rWiNew = math::clamp(rWiNew, 0, m_wiResolution - 1);

			for (int dcWi = 0; dcWi < 2; dcWi++) {
				int cWiNew = cWi + dcWi;
// 				if (cWiNew < 0 || cWiNew >= m_wiResolution)
// 					continue;
				cWiNew = math::clamp(cWiNew, 0, m_wiResolution - 1);
				
				Vector tableWi;
				tableWi.x = (cWiNew + 0.5) / m_wiResolution * 2.0 - 1.0;
				tableWi.y = ((m_wiResolution - rWiNew - 1) + 0.5) / m_wiResolution * 2.0 - 1.0;
				double tmp = tableWi.x * tableWi.x + tableWi.y * tableWi.y;
				if (tmp > 1.0)
					continue;
				tableWi.z = std::sqrt(1.0 - tmp);

				for (int drWo = 0; drWo < 2; drWo++) {
					int rWoNew = rWo + drWo;
// 					if (rWoNew < 0 || rWoNew >= m_woResolution)
// 						continue;
					rWoNew = math::clamp(rWoNew, 0, m_woResolution - 1);

					for (int dcWo = 0; dcWo < 2; dcWo++) {
						int cWoNew = cWo + dcWo;
// 						if (cWoNew < 0 || cWoNew >= m_woResolution)
// 							continue;
						cWoNew = math::clamp(cWoNew, 0, m_woResolution - 1);

						Vector tableWo;
						tableWo.x = (cWoNew + 0.5) / m_woResolution * 2.0 - 1.0;
						tableWo.y = ((m_woResolution - rWoNew - 1) + 0.5) / m_woResolution * 2.0 - 1.0;
						tmp = tableWo.x * tableWo.x + tableWo.y * tableWo.y;
						if (tmp > 1.0)
							continue;
						tableWo.z = std::sqrt(1.0 - tmp);

						Spectrum tmpValue = m_angularScales->getPixel(Point2i(rWoNew * m_woResolution + cWoNew,
							rWiNew * m_wiResolution + cWiNew));
						double dCosine = std::max(0.0, dot(wiMacro, tableWi)) * std::max(0.0, dot(woMacro, tableWo));
						double tmpWeight = math::fastexp(1000.0 * (dCosine - 1.0));

						res += tmpValue * tmpWeight;
						weights += tmpWeight;
					}
				}
			}
		}

		if (weights < 1e-5) {
// 			Log(EInfo, "wi = [%d, %d](%.6f, %.6f, %.6f), wo = [%d, %d](%.6f, %.6f, %.6f), (%.6f, %.6f, %.6f)", 
// 				rWi, cWi, 
// 				wiMacro.x, wiMacro.y, wiMacro.z,
// 				rWo, cWo,
// 				woMacro.x, woMacro.y, woMacro.z,
// 				res[0], res[1], res[2]);

			res = Spectrum(1.0f);
		} else {
			res /= weights;
		}

// 		Log(EInfo, "wi = (%.6f, %.6f, %.6f), wo = (%.6f, %.6f, %.6f), %.6f", 
// 			wiMacro.x, wiMacro.y, wiMacro.z,
// 			woMacro.x, woMacro.y, woMacro.z, res[0]);
		return res;
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		Spectrum s(1.f), t(1.f);
		if (m_angularScales != NULL)
			s = evalScale(bRec);
		if (m_spatialScales != NULL)
			t = evalSpatialScale(bRec);
		if (s.isZero() || t.isZero())
			return Spectrum(0.f);
		Spectrum spec = m_bsdf->eval(bRec, measure);
		return spec * s * t;
	}

	inline int getScaleMatrixRow(const BSDFSamplingRecord &bRec) const {
		Vector wiWorld = bRec.its.toWorld(bRec.wi);
		Vector wiMacro = bRec.its.baseFrame.toLocal(wiWorld);
		Point2 wiTex = warp::uniformHemisphereToSquareConcentric(wiMacro);
		
		// piecewise constant
		int wiNumCells = m_wiResolution;
		int c1 = math::clamp(math::floorToInt(wiTex.x * wiNumCells), 0, wiNumCells - 1);
		int r1 = math::clamp(math::floorToInt(wiTex.y * wiNumCells), 0, wiNumCells - 1);
		return r1 * m_wiResolution + c1;
	}

	inline int getScaleMatrixCol(const BSDFSamplingRecord &bRec) const {
		Vector woWorld = bRec.its.toWorld(bRec.wo);
		Vector woMacro = bRec.its.baseFrame.toLocal(woWorld);
		Point2 woTex = warp::uniformHemisphereToSquareConcentric(woMacro);

		// piecewise constant
		int woNumCells = m_woResolution;
		int c2 = math::clamp(math::floorToInt(woTex.x * woNumCells), 0, woNumCells - 1);
		int r2 = math::clamp(math::floorToInt(woTex.y * woNumCells), 0, woNumCells - 1);
		return r2 * m_woResolution + c2;
	}

	inline Spectrum sampleScaleMatrixCol(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		Point2i p;
		p.y = getScaleMatrixRow(bRec);
		p.x = m_rowCondCdfs[p.y].sample(sample.x);
		pdf = m_rowCondCdfs[p.y][p.x];
		Spectrum spec = m_angularScales->getPixel(p);

		// reconstruct wo
		ref<Sampler> sampler = m_samplers[Thread::getID() % 233];
		int c2 = p.x % m_woResolution;
		int r2 = p.x / m_woResolution;
		Point2 woTex((c2 + sampler->next1D()) / (double)m_woResolution,
			(r2 + sampler->next1D()) / (double)m_woResolution);
		Vector woMacro = warp::squareToUniformHemisphereConcentric(woTex);
		Vector woWorld = bRec.its.baseFrame.toWorld(woMacro);
		
		bRec.wo = bRec.its.toLocal(woWorld);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		return spec;
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (m_useMIS) {
			Float pdf1 = m_bsdf->pdf(bRec, measure);
			Point2i p;
			p.x = getScaleMatrixCol(bRec);
			p.y = getScaleMatrixRow(bRec);
			Float pdf2 = m_rowCondCdfs[p.y][p.x];
			return 0.5 * (pdf1 + pdf2);
		} else {
			return m_bsdf->pdf(bRec, measure);
		}
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (m_useMIS) {
			ref<Sampler> sampler = m_samplers[Thread::getID() % 233];
			if (sampler->next1D() <= 0.5) {
				Spectrum spec = m_bsdf->sample(bRec, pdf, sample);
				if (spec.isZero())
					return Spectrum(0.f);
				Spectrum scale = evalScale(bRec);

				Point2i p;
				p.y = getScaleMatrixRow(bRec);
				p.x = getScaleMatrixCol(bRec);
				Float pdf2 = m_rowCondCdfs[p.y][p.x];
				Float weight = miWeight(pdf, pdf2);
				return spec * scale * weight * 2.0;
			} else {
				Spectrum scale = sampleScaleMatrixCol(bRec, pdf, sample);
				Spectrum spec = m_bsdf->eval(bRec);
				if (spec.isZero())
					return Spectrum(0.f);

				Float pdf2 = m_bsdf->pdf(bRec);
				Float weight = miWeight(pdf, pdf2);
				return spec * scale * weight * 2.0;
			}
		} else {
			// naive cosine sampling
// 			bRec.wo = warp::squareToCosineHemisphere(sample);
// 			bRec.sampledComponent = 0;
// 			bRec.sampledType = EDiffuseReflection;
// 			pdf = warp::squareToCosineHemispherePdf(bRec.wo);
// 			return eval(bRec, ESolidAngle) / pdf;

			Spectrum spec = m_bsdf->sample(bRec, pdf, sample);
			if (spec.isZero())
				return Spectrum(0.f);
			
			Spectrum s(1.f), t(1.f);
			
			if (m_angularScales != NULL) {
				if (!m_useGrad)
					s = evalScale(bRec);
				else
					s = evalScaleAndGrad(bRec);
			}

			if (m_spatialScales != NULL) {
				t = evalSpatialScale(bRec);
			}

// 			if (!spec.isZero() && !s.isZero()) {
// 				Vector wiWorld = bRec.its.toWorld(bRec.wi);
// 				Vector wiMacro = bRec.its.baseFrame.toLocal(wiWorld);
// 				Vector woWorld = bRec.its.toWorld(bRec.wo);
// 				Vector woMacro = bRec.its.baseFrame.toLocal(woWorld);
// 				Log(EInfo, "======================");
// 				Log(EInfo, "wiMacro = (%.6f, %.6f, %.6f)", wiMacro.x, wiMacro.y, wiMacro.z);
// 				Log(EInfo, "woMacro = (%.6f, %.6f, %.6f)", woMacro.x, woMacro.y, woMacro.z);
// 				Log(EInfo, "bsdfWeight = (%.6f, %.6f, %.6f)", spec[0], spec[1], spec[2]);
// 				Log(EInfo, "S = (%.6f, %.6f, %.6f)", s[0], s[1], s[2]);
// 				Log(EInfo, "pdf = %.6f", pdf);
// 			}

			return spec * s * t;
		}
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sp) const {
		Float pdf;
		return sample(bRec, pdf, sp);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
			m_bsdf = static_cast<BSDF *>(child);
		}
		else {
			BSDF::addChild(name, child);
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "TabulatedScaledBSDF[" << endl
			<< "  filename = \"" << m_angularScaleFilename << "\"," << endl;
		oss << "]";
		return oss.str();
	}

	inline Point2 transformUV(const Point2 &_uv) const {
		Point2 uv(_uv);
		uv.x *= m_uvScale.x;
		uv.y *= m_uvScale.y;
		uv.x = uv.x - math::floorToInt(uv.x);
		uv.y = uv.y - math::floorToInt(uv.y);
		return uv;
	}

	Float getRoughness(const Intersection &its, int component) const {
		return m_bsdf->getRoughness(its, component);
	}

	Spectrum getLobeAlbedo(const Intersection &its, int component) const {
		return m_bsdf->getLobeAlbedo(its, component);
	}

	Float getLobeRoughness(const Intersection &its, int component) const {
		return m_bsdf->getLobeRoughness(its, component);
	}

	inline Float miWeight(Float pdfA, Float pdfB) const {
		pdfA *= pdfA;
		pdfB *= pdfB;
		return pdfA / (pdfA + pdfB);
	}

	MTS_DECLARE_CLASS()
public:
	ref<Bitmap> m_angularScales;
	ref<Bitmap> m_spatialScales;
	ref<BSDF> m_bsdf;
	std::string m_angularScaleFilename;
	std::string m_spatialScaleFilename;
	std::string m_spatialInterp;
	Vector2i m_lobeSize;
	
	Vector2i m_xyResolution;
	Vector2 m_uvScale;

	int m_wiResolution;
	int m_woResolution;
	
	bool m_wiUseFullSphere;
	bool m_woUseFullSphere;
	bool m_useGrad;

	bool m_useMIS;
	ref_vector<Sampler> m_samplers;
	std::vector<DiscreteDistribution> m_rowCondCdfs;
};

MTS_IMPLEMENT_CLASS_S(TabulatedScaledBSDF, false, BSDF)
MTS_EXPORT_PLUGIN(TabulatedScaledBSDF, "Tabulated scaled BSDF");
MTS_NAMESPACE_END
