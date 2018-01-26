#if !defined(__MICROFACET_MULTI_H)
#define __MICROFACET_MULTI_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/properties.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

static inline double abgam(double x) {
	double gam[10], temp;
	gam[0] = 1./ 12.;
	gam[1] = 1./ 30.;
	gam[2] = 53./ 210.;
	gam[3] = 195./ 371.;
	gam[4] = 22999./ 22737.;
	gam[5] = 29944523./ 19733142.;
	gam[6] = 109535241009./ 48264275462.;
	temp = 0.5*log (2*M_PI) - x + (x - 0.5)*log (x)
		+ gam[0]/(x + gam[1]/(x + gam[2]/(x + gam[3]/(x + gam[4] /
		(x + gam[5]/(x + gam[6]/x))))));

	return temp;
}

static inline double gamma(double x) {
	double result;
	result = math::fastexp(abgam(x + 5))/(x*(x + 1)*(x + 2)*(x + 3)*(x + 4));
	return result;
}

static inline double beta(double m, double n) {
	return (gamma(m)*gamma(n)/gamma(m + n));
}

#define vec3 Vector
#define vec2 Vector

struct RayInfo {
	// direction
	vec3 w;
	vec3 mesoN;
	Float theta;
	Float cosTheta;
	Float sinTheta;
	Float tanTheta;
	Float alpha;
	Float Lambda;

	void updateDirection(const vec3& w, const Float alpha_x, const Float alpha_y) {
		this->w = w;
		mesoN = Normal(0.0, 0.0, 1.0);
		theta = acos(w.z);
		cosTheta = w.z;
		sinTheta = sin(theta);
		tanTheta = sinTheta / cosTheta;
		const Float invSinTheta2 = 1.0f / (1.0f - w.z*w.z);
		const Float cosPhi2 = w.x*w.x*invSinTheta2;
		const Float sinPhi2 = w.y*w.y*invSinTheta2;
		alpha = sqrt(cosPhi2*alpha_x*alpha_x + sinPhi2*alpha_y*alpha_y); 
		// Lambda
		if(w.z > 0.9999f)
			Lambda = 0.0f;
		else if(w.z < -0.9999f)
			Lambda = -1.0f;
		else {	
			const Float a = 1.0f/tanTheta/alpha;
			Lambda = 0.5f*(-1.0f + ((a>0)?1.0f:-1.0f) * sqrt(1 + 1/(a*a)));
		}
	}

	void updateDirection(const Vector &w, const Spectrum &moments0, 
		Float sigmaX2, Float sigmaY2, Float cxy) {
		this->w = w;
		mesoN = normalize(Normal(-moments0[0], -moments0[1], 1.0));
		theta = acos(w.z);
		cosTheta = w.z;
		sinTheta = sin(theta);
		tanTheta = sinTheta / cosTheta;
		const Float invSinTheta2 = 1.0 / (1.0 - w.z * w.z);
		const Float cosPhi2 = w.x * w.x * invSinTheta2;
		const Float sinPhi2 = w.y * w.y * invSinTheta2;
		const Float cosPhi = sqrt(cosPhi2);
		const Float sinPhi = sqrt(sinPhi2);
		alpha = cosPhi2 * sigmaX2 + sinPhi2 * sigmaY2 + 2.0 * cosPhi * sinPhi * cxy;
		alpha = sqrt(2.0 * alpha);
		if (w.z > 0.9999f)
			Lambda = 0.0f;
		else if (w.z < -0.9999f)
			Lambda = -1.0f;
		else {
			const Float muPhi = cosPhi * moments0[0] + sinPhi * moments0[1];
			const Float a = (1.0 / tanTheta - muPhi) / alpha;
			Lambda = 0.5 * (math::erf(a) - 1.0f) + math::fastexp(-a * a) / a * 0.5 * sqrt(INV_PI);
		}
	}

	// height
	Float h;
	Float C1;
	Float G1;

	void updateHeight(const Float& h) {
		this->h = h;
		C1 = std::min(1.0, std::max(0.0, 0.5*(h+1.0)));

		if(this->w.z > 0.9999f)
			G1 = 1.0f;
		else if(this->w.z <= 0.0f)
			G1 = 0.0f;
		else
			G1 = pow(this->C1, this->Lambda / this->mesoN.z);
	}
};

inline Float invC1(const Float U) {
	const Float h = std::max(-1.0, std::min(1.0, 2.0*U-1.0));
	return h;	
}

inline Float sampleHeight(const RayInfo& ray, const Float U) {
	if(ray.w.z > 0.9999f)
		return std::numeric_limits<Float>::max();
	if(ray.w.z < -0.9999f) {
		const Float value = invC1(U*ray.C1);
		return value;
	}
	if(fabs(ray.w.z) < 0.0001f)
		return ray.h;

	// probability of intersection
	if (U > 1.0f - ray.G1) // leave the microsurface
		return std::numeric_limits<Float>::max();

	const Float h = invC1( 
		ray.C1 / pow((1.0f-U),ray.mesoN.z/ray.Lambda)
		);
	return h;
}

Float D_ggx(const vec3& wm, const Float alpha_x, const Float alpha_y) {
	if( wm.z <= 0.0f)
		return 0.0f;

	// slope of wm
	const Float slope_x = -wm.x/wm.z;
	const Float slope_y = -wm.y/wm.z;

	// P22
	const Float tmp = 1.0f + slope_x*slope_x/(alpha_x*alpha_x) + slope_y*slope_y/(alpha_y*alpha_y);
	const Float P22 = 1.0f / (M_PI * alpha_x * alpha_y) / (tmp * tmp);

	// value
	const Float value = P22 / (wm.z*wm.z*wm.z*wm.z);
	return value;
}

vec2 sampleP22_11(const Float theta_i, const Float U, const Float U_2, const Float alpha_x, const Float alpha_y) {
	vec2 slope;

	if(theta_i < 0.0001f) {
		const Float r = sqrt(U/(1.0f-U));
		const Float phi = 6.28318530718f * U_2;
		slope.x = r * cos(phi);
		slope.y = r * sin(phi);
		return slope;
	}

	// constant
	const Float sin_theta_i = sin(theta_i);
	const Float cos_theta_i = cos(theta_i);
	const Float tan_theta_i = sin_theta_i/cos_theta_i;

	// slope associated to theta_i
	const Float slope_i = cos_theta_i/sin_theta_i;

	// projected area
	const Float projectedarea = 0.5f * (cos_theta_i + 1.0f);
	if(projectedarea < 0.0001f || projectedarea!=projectedarea)
		return vec2(0,0,0);
	// normalization coefficient
	const Float c = 1.0f / projectedarea;

	const Float A = 2.0f*U/cos_theta_i/c - 1.0f;
	const Float B = tan_theta_i;
	const Float tmp = 1.0f / (A*A-1.0f);

	const Float D = sqrt(std::max(0.0, B*B*tmp*tmp - (A*A-B*B)*tmp));
	const Float slope_x_1 = B*tmp - D;
	const Float slope_x_2 = B*tmp + D;
	slope.x = (A < 0.0f || slope_x_2 > 1.0f/tan_theta_i) ? slope_x_1 : slope_x_2;

	Float U2;
	Float S;
	if(U_2 > 0.5f) {
		S = 1.0f;
		U2 = 2.0f*(U_2-0.5f);
	}
	else {
		S = -1.0f;
		U2 = 2.0f*(0.5f-U_2);
	}
	const Float z = (U2*(U2*(U2*0.27385f-0.73369f)+0.46341f)) / (U2*(U2*(U2*0.093073f+0.309420f)-1.000000f)+0.597999f);
	slope.y = S * z * sqrt(1.0f+slope.x*slope.x);

	return slope;
}

vec3 sampleVNDF(const vec3& wi, const Float alpha_x, const Float alpha_y, Sampler *sampler) {
	const Float U1 = sampler->next1D();
	const Float U2 = sampler->next1D();

	// sample D_wi

	// stretch to match configuration with alpha=1.0	
	const vec3 wi_11 = normalize(vec3(alpha_x * wi.x, alpha_y * wi.y, wi.z));

	// sample visible slope with alpha=1.0
	vec2 slope_11 = sampleP22_11(acos(wi_11.z), U1, U2, alpha_x, alpha_y);

	// align with view direction
	const Float phi = atan2(wi_11.y, wi_11.x);
	vec2 slope(cos(phi)*slope_11.x - sin(phi)*slope_11.y, sin(phi)*slope_11.x + cos(phi)*slope_11.y, 0);

	// stretch back
	slope.x *= alpha_x;
	slope.y *= alpha_y;

	// if numerical instability
	if ((slope.x != slope.x) || !std::isfinite(slope.x)) {
		if (wi.z > 0) return vec3(0.0f,0.0f,1.0f);
		else return normalize(vec3(wi.x, wi.y, 0.0f));
	}

	// compute normal
	const vec3 wm = normalize(vec3(-slope.x, -slope.y, 1.0f));

	return wm;
}

MTS_NAMESPACE_END

#endif
