#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/mipmap.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

class SHScaledBSDF : public BSDF {
public:
	SHScaledBSDF(const Properties &props) : BSDF(props) {
		// avoid negative value
		m_coeffsOffset = 100.0;

		m_numCoeffs = props.getInteger("numCoeffs", 1);

		Float uvscale = props.getFloat("uvscale", 1.0f);
		m_uvScale = Vector2(
			props.getFloat("uscale", uvscale),
			props.getFloat("vscale", uvscale)
		);

		m_coeffsFilenamePrefix = props.getString("prefix", "");
	}

	SHScaledBSDF(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		m_numCoeffs = stream->readInt();
		m_uvScale = Vector2(stream->readFloat(), stream->readFloat());
		m_coeffsFilenamePrefix = stream->readString();

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeInt(m_numCoeffs);
		stream->writeFloat(m_uvScale.x);
		stream->writeFloat(m_uvScale.y);
		stream->writeString(m_coeffsFilenamePrefix);
	}

	void configure() {
		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide | ESpatiallyVarying);

		m_usesRayDifferentials = false;

		// load textures
		m_shCoeffs.resize(m_numCoeffs);
		for (int i = 0; i < m_numCoeffs; i++) {
			Properties props("bitmap");
			props.setString("wrapMode", "repeat");
			props.setFloat("gamma", 1.0f);
			props.setFloat("uscale", m_uvScale.x);
			props.setFloat("vscale", m_uvScale.y);
			props.setString("filterType", "nearest");
			std::ostringstream oss;
			oss << m_coeffsFilenamePrefix << i << ".exr";
			props.setString("filename", oss.str());
			m_shCoeffs[i] = static_cast<Texture *>(PluginManager::getInstance()->
				createObject(MTS_CLASS(Texture), props));
			m_shCoeffs[i]->configure();
		}
		Log(EInfo, "Load %d SH coeffcients as textures", m_numCoeffs);

		BSDF::configure();

		/*
		Intersection its;
		its.uv.x = 6; its.uv.y = 2;
		its.uv.x = its.uv.x / 16.0 / 40.0;
		its.uv.y = its.uv.y / 16.0 / 40.0;
		Spectrum spec = m_shCoeffs[0]->eval(its, false);
		Log(EInfo, "%.6f, %.6f, %.6f", spec[0], spec[1], spec[2]);
		*/

		/*
		SHVector shVec(3);
		shVec(0, 0) = 2.5429;
		shVec(1, -1) = 1.4690;
		shVec(1, 0) = -2.1013;
		shVec(1, 1) = -1.4717;
		shVec(2, -2) = 0.7166;
		shVec(2, -1) = -0.5384;
		shVec(2, 0) = 2.3998;
		shVec(2, 1) = 0.8145;
		shVec(2, 2) = 0.0227;

		int m_size = 64;
		std::vector<std::vector<Vector3d> > bsdfValues(m_size, std::vector<Vector3d>(m_size));
		for (int i = 0; i < m_size; i++)
			for (int j = 0; j < m_size; j++)
				bsdfValues[i][j] = Vector3d(0.0f);

		for (int i = 0; i < m_size; i++) {
			for (int j = 0; j < m_size; j++) {
				Vector3d res(0.0);
				double dx = (j + 0.5) / (double)m_size * 2.0 - 1.0;
				double dy = (i + 0.5) / (double)m_size * 2.0 - 1.0;
				if (dx * dx + dy * dy > 1.0)
					continue;
				double dz = sqrt(1.0 - dx * dx - dy * dy);

				for (int c = 0; c < 3; c++) {
					res[c] += shVec.evalXYZ(Vector(dx, dy, dz));
				}

				res /= 1.0;

				for (int c = 0; c < 3; c++) {
					bsdfValues[m_size - i - 1][j][c] = res[c];
				}
			}
		}

		ref<Bitmap> bitmap = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, Vector2i(m_size));
		float *data = bitmap->getFloat32Data();
		for (int r = 0; r < m_size; r++) {
			for (int c = 0; c < m_size; c++) {
				*data++ = bsdfValues[r][c][0];
				*data++ = bsdfValues[r][c][1];
				*data++ = bsdfValues[r][c][2];
			}
		}
		ref<FileStream> stream = new FileStream("A_sh_9_mitsuba.exr", FileStream::ETruncWrite);
		bitmap->write(Bitmap::EOpenEXR, stream);
		*/
	}

	Spectrum evalScale(const BSDFSamplingRecord &bRec) const {
		int bands = (int)std::sqrt((Float)m_numCoeffs + 1e-4);
		Spectrum res(0.f);
		std::vector<SHVector> shVecs(3, SHVector(bands));

		int idx = 0;
		for (int l = 0; l < bands; l++) {
			for (int m = -l; m <= l; m++) {
				Spectrum coeff = m_shCoeffs[idx]->eval(bRec.its, false);
				for (int c = 0; c < SPECTRUM_SAMPLES; c++) {
					shVecs[c](l, m) = coeff[c] - m_coeffsOffset;
				}
				idx++;
			}
		}

		Vector woWorld = bRec.its.toWorld(bRec.wo);
		Vector woMacro = bRec.its.baseFrame.toLocal(woWorld);
		for (int c = 0; c < SPECTRUM_SAMPLES; c++) {
			res[c] = shVecs[c].evalXYZ(woMacro);
			res[c] = std::max(0.0, std::min(1.0, res[c]));
		}

		/*
		if (woMacro.z > 0.5) {
		Log(EInfo, "=================");
		for (int l = 0; l < bands; l++) {
			for (int m = -l; m <= l; m++) {
				Log(EInfo, "%.6f", shVecs[0](l, m));
			}
		}
		Log(EInfo, "SHeval = %.6f", shVecs[0].evalXYZ(woMacro));
		Log(EInfo, "v = (%.6f, %.6f, %.6f), %.6f, %.6f, %.6f", 
			woMacro.x, woMacro.y, woMacro.z,
			res[0], res[1], res[2]);
		}
		*/

		return res;
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		Spectrum scale = evalScale(bRec);
		if (scale.isZero())
			return Spectrum(0.f);
		Spectrum spec = m_bsdf->eval(bRec, measure);
		return spec * scale;
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		return m_bsdf->pdf(bRec, measure);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		Spectrum spec = m_bsdf->sample(bRec, pdf, sample);
		if (spec.isZero())
			return Spectrum(0.f);
		Spectrum scale = evalScale(bRec);
		return spec * scale;
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
		oss << "SHScaledBSDF[" << endl
			<< "  filename = \"" << m_coeffsFilenamePrefix << "\"," << endl;
		oss << "]";
		return oss.str();
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

	MTS_DECLARE_CLASS()
public:
	int m_numCoeffs;
	Float m_coeffsOffset;
	std::string m_coeffsFilenamePrefix;
	Vector2 m_uvScale;
	ref_vector<Texture> m_shCoeffs;
	ref<BSDF> m_bsdf;
};

MTS_IMPLEMENT_CLASS_S(SHScaledBSDF, false, BSDF)
MTS_EXPORT_PLUGIN(SHScaledBSDF, "SH scaled BSDF");
MTS_NAMESPACE_END
