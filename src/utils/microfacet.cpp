#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include "../bsdfs/microfacet.h"

MTS_NAMESPACE_BEGIN

class TestMicrofacet : public Utility {
	int run(int argc, char **argv) {
		MicrofacetDistribution dist(MicrofacetDistribution::EBeckmann, 1.0f);
		Vector m(0, 0, 1);
		
		Vector wi(std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3]));
		wi = normalize(wi);

		int nTheta = 180;
		int nPhi = 360;
		double dTheta = 0.5 * M_PI / (double)nTheta;
		double dPhi = 2.0 * M_PI / (double)nPhi;

		double totVisible = 0.0;
		for (int i = 0; i < nTheta; i++) {
			double theta = (i + 0.5f) * dTheta;
			double sinTheta = std::sin(theta);
			for (int j = 0; j < nPhi; j++) {
				double phi = (j + 0.5f) * dPhi;
				Vector wm(sinTheta * std::cos(phi), sinTheta * std::sin(phi), std::cos(theta));
				double cosTerm = std::max(0.0f, dot(wi, wm));
				totVisible += dist.eval(wm) * cosTerm * sinTheta * dTheta * dPhi;
			}
		}
		Log(EInfo, "total visible = %.6f", totVisible);

		double res = 0.0;
		for (int i = 0; i < nTheta; i++) {
			double theta = (i + 0.5f) * dTheta;
			double sinTheta = std::sin(theta);
			for (int j = 0; j < nPhi; j++) {
				double phi = (j + 0.5f) * dPhi;
				Vector wo(sinTheta * std::cos(phi), sinTheta * std::sin(phi), std::cos(theta));
				Vector wh = wi + wo;
				wh = normalize(wh);
				double cosTerm = std::max(0.0f, dot(wi, wh));
				double visD_wh = dist.eval(wh) * cosTerm / totVisible;
				double shadowGivenMaskTerm = (1.0 + dist.smithLambda(wi, m)) / (1.0 + dist.smithLambda(wi, m) + dist.smithLambda(wo, m));
				//double shadowGivenMaskTerm = dist.smithG1(wo, m);
				res += shadowGivenMaskTerm * visD_wh * sinTheta * dTheta * dPhi / (4.0 * wh.z);
			}
		}
		Log(EInfo, "G1 = %.6f", dist.smithG1(wi, m));
		Log(EInfo, "avg G2 = %.6f %.6f", res, res * dist.smithG1(wi, m));
		return 0;
	}

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(TestMicrofacet, "Play with microfacet distribution")
MTS_NAMESPACE_END
