#ifndef FILTER_INCLUDED
#define FILTER_INCLUDED
#define PI_2 6.28318530718
#define SIGMA_SQUARED 2
#include <algorithm>

class ContinuousFilter{
public:
	virtual double value(double u, double v) = 0;
	int support_h;
	int support_w;
};

class DerivativeFilter : public ContinuousFilter{
public:
	virtual double value(double u, double v) = 0;
	virtual double value_dx(double u, double v) = 0;
	virtual double value_dy(double u, double v) = 0;
	virtual double value_dxdy(double u, double v) = 0;
	virtual double value_dxdx(double u, double v) = 0;
	virtual double value_dydy(double u, double v) = 0;
	virtual double value_dxdx_dydy(double u, double v) = 0;
	int support_h;
	int support_w;
};

class BilinearFilter : public ContinuousFilter{
public:
	BilinearFilter(){
		support_h = 2;
		support_w = 2;
	}
	double value(double u, double v){
		return std::max<double>(1.f - abs(u), 0.f)*std::max<double>(1.f - abs(v), 0.f);
	}
};

class DiscreteFilter{
public:
	virtual double value(int u, int v) = 0;
	int radius_h;
	int radius_w;
};

class Identity_Filter : public DiscreteFilter{
public:
	Identity_Filter(){
		radius_h = 0;
		radius_w = 0;
	}
	double  value(int u, int v){
		if (v == 0 && u == 0){
			return 1.f;
		}
		return 0.f;
	}
};

class Derivative_w : public DiscreteFilter{
public:
	Derivative_w(){
		radius_h = 0;
		radius_w = 1;
	}
	double  value(int u, int v){
		if (v == 0){
			if (u == -1) return -1.f;
			if (u == 0)  return 1.f;
		}
		return 0.f;
	}
};

class Derivative_ww : public DiscreteFilter{
public:
	Derivative_ww(){
		radius_h = 0;
		radius_w = 1;
	}
	double  value(int u, int v){
		if (v == 0){
			if (u == -1) return 1.f;
			if (u == 0)  return -2.f;
			if (u == 1)  return 1.f;
		}
		return 0.f;
	}
};

class Derivative_h : public DiscreteFilter{
public:
	Derivative_h(){
		radius_h = 1;
		radius_w = 0;
	}
	double  value(int u, int v){
		if (u == 0){
			if (v == -1) return -1.f;
			if (v == 0)  return 1.f;
		}
		return 0.f;
	}
};

class Derivative_hh : public DiscreteFilter{
public:
	Derivative_hh(){
		radius_h = 1;
		radius_w = 0;
	}
	double  value(int u, int v){
		if (u == 0){
			if (v == -1) return 1.f;
			if (v == 0)  return -2.f;
			if (v == 1)  return 1.f;
		}
		return 0.f;
	}
};

class Derivative_hw : public DiscreteFilter{
public:
	Derivative_hw(){
		radius_h = 1;
		radius_w = 1;
	}
	double  value(int u, int v){
		if (u == -1){
			if (v == -1) return 1.f;
			if (v == 0)  return -1.f;
		}
		if (u == 0){
			if (v == -1) return -1.f;
			if (v == 0)  return 1.f;
		}
		return 0.f;
	}
};

class Laplacian : public DiscreteFilter{
public:
	Laplacian(){
		radius_h = 1;
		radius_w = 1;
	}
	double value(int u, int v){
		if (u == -1){
			if (v == 0)  return 1.f;
		}
		if (u == 0){
			if (v == -1) return 1.f;
			if (v == 0)  return -4.f;
			if (v == 1) return 1.f;
		}
		if (u == 1){
			if (v == 0)  return 1.f;
		}
		return 0.f;
	}
};

class ContinuousBspline3 : public ContinuousFilter {
public:
	ContinuousBspline3(){
		support_h = 4;
		support_w = 4;
	}
	double value(double r)
	{
		r = abs(r);
		if (r < 1.f) return (4.f + r*r*(-6.f + 3.f*r)) / 6.f;
		else if (r < 2.f) return  (8.f + r*(-12.f + (6.f - r)*r)) / 6.f;
		else return 0.f;

	}
	double value(double u, double v){
		return value(u)*value(v);
	}
};


class ContinuousBspline3_Derivative : public DerivativeFilter {
public:
	ContinuousBspline3_Derivative(){
		support_h = 4;
		support_w = 4;
	}

	double value(double r)
	{
		r = abs(r);
		if (r < 1.f) return (4.f + r*r*(-6.f + 3.f*r)) / 6.f;
		else if (r < 2.f) return  (8.f + r*(-12.f + (6.f - r)*r)) / 6.f;
		else return 0.f;
	}

	double value_d(double r)
	{
		double sign_r = r > 0.f ? 1.f : -1.f;
		r = abs(r);
		if (r < 1.f){
			return (r*(r*3.f - 4.f) / 2.f)*sign_r;
		}
		else if (r < 2.f){
			return ((r*(4.f - r) - 4.f) / 2.f)*sign_r;
		}
		else return 0.f;
	}

	double value_dd(double r)
	{
		r = abs(r);
		if (r < 1.f){
			return (3.f*r-2.f);
		}
		else if (r < 2.f){
			return 2.f -r;
		}
		else return 0.f;
	}

	double value(double u, double v){
		return value(u)*value(v);
	}
	double value_dx(double u, double v){
		return value_d(u)*value(v);
	}
	double value_dy(double u, double v){
		return value(u)*value_d(v);
	}
	double value_dxdy(double u, double v){
		return value_d(u)*value_d(v);
	}
	double value_dxdx(double u, double v){
		return value_dd(u)*value(v);
	}
	double value_dydy(double u, double v){
		return value(u)*value_dd(v);
	}
	double value_dxdx_dydy(double u, double v){
		return value_dd(u)*value(v) + value(u)*value_dd(v);
	}
};


class Gaussian2_Derivatives{
public:
	static double value(double u, double v){
		return exp(-(u*u + v*v) / (2.f*SIGMA_SQUARED)) / (PI_2*SIGMA_SQUARED);
	}
	static double value_dx(double u, double v){
		return -(u / SIGMA_SQUARED)*value(u,v);
	}
	static double value_dy(double u, double v){
		return -(v / SIGMA_SQUARED)*value(u, v);
	}
	static double value_dxdy(double u, double v){
		return ((u*v)/(SIGMA_SQUARED*SIGMA_SQUARED))*value(u, v);
	}
	static double value_dxdx(double u, double v){
		double val = value(u, v);
		return ((u*u - SIGMA_SQUARED)/(SIGMA_SQUARED*SIGMA_SQUARED))*val;
	}
	static double value_dydy(double u, double v){
		double val = value(u, v);
		return ((v*v - SIGMA_SQUARED) / (SIGMA_SQUARED*SIGMA_SQUARED))*val;
	}
	static double value_dxdx_dydy(double u, double v){
		double val = value(u, v);
		return ((u*u + v*v - 2.f*SIGMA_SQUARED) / (SIGMA_SQUARED*SIGMA_SQUARED))*val;
	}
};

class Bspline3_Derivatives{
public:
	static double value(double r)
	{
		r = abs(r);
		if (r < 1.f) return (4.f + r*r*(-6.f + 3.f*r)) / 6.f;
		else if (r < 2.f) return  (8.f + r*(-12.f + (6.f - r)*r)) / 6.f;
		else return 0.f;
	}

	static double value_d(double r)
	{
		double sign_r = r > 0.f ? 1.f : -1.f;
		r = abs(r);
		if (r < 1.f){
			return (r*(r*3.f - 4.f) / 2.f)*sign_r;
		}
		else if (r < 2.f){
			return ((r*(4.f - r) - 4.f) / 2.f)*sign_r;
		}
		else return 0.f;
	}

	static double value_dd(double r)
	{
		r = abs(r);
		if (r < 1.f){
			return (3.f*r - 2.f);
		}
		else if (r < 2.f){
			return 2.f - r;
		}
		else return 0.f;
	}

	static double value(double u, double v){
		return value(u)*value(v);
	}
	static double value_dx(double u, double v){
		return value_d(u)*value(v);
	}
	static double value_dy(double u, double v){
		return value(u)*value_d(v);
	}
	static double value_dxdy(double u, double v){
		return value_d(u)*value_d(v);
	}
	static double value_dxdx(double u, double v){
		return value_dd(u)*value(v);
	}
	static double value_dydy(double u, double v){
		return value(u)*value_dd(v);
	}
	static double value_dxdx_dydy(double u, double v){
		return value_dd(u)*value(v) + value(u)*value_dd(v);
	}
};

class Bspline5_Derivatives{
public:
	static double value(double r)
	{
		r = abs(r);
		if (r < 1.f) return (66.f + r*r*(-60.f
			+ (30.f - 10.f*r)*r*r)) / 120.f;
		else if (r < 2.f) return (51.f + r*(75.f + r*(-210.f
			+ r*(150.f + r*(-45.f + 5.f*r))))) / 120.f;
		else if (r < 3.f) return (243.f + r*(-405.f + r*(270.f
			+ r*(-90.f + (15.f - r)*r)))) / 120.f;
		else return 0.f;
	}
	static double value_d(double r)
	{
		double sign_r = r > 0.f ? 1.f : -1.f;
		r = abs(r);
		if (r < 1.f){
			return  (r*(-1.f + (r*r*(1.f - (5.f*r/ 12.f)))))*sign_r;
		}
		else if (r < 2.f){
			return  ((15.f + r*(-84.f + r*(90.f + r*(-36.f + 5.f*r))))/24.f)*sign_r;
		}
		else if (r < 3.f) return (-(r - 3.f)*(r - 3.f)*(r - 3.f)*(r - 3.f)/ 24.f)*sign_r;
		else return 0.f;
	}

	static double value_dd(double r)
	{
		r = abs(r);
		if (r < 1.f){
			return  -1.f + r*r*(3.f  - (5.f*r/ 3.f));
		}
		else if (r < 2.f){
			return  (-21.f + r*(45 + r*(-27.f + 5*r)))/6.f;
		}
		else if (r < 3.f) return  (-(r - 3.f)*(r - 3.f)*(r - 3.f) / 6.f);
		else return 0.f;
	}

	static double value(double u, double v){
		return value(u)*value(v);
	}
	static double value_dx(double u, double v){
		return value_d(u)*value(v);
	}
	static double value_dy(double u, double v){
		return value(u)*value_d(v);
	}
	static double value_dxdy(double u, double v){
		return value_d(u)*value_d(v);
	}
	static double value_dxdx(double u, double v){
		return value_dd(u)*value(v);
	}
	static double value_dydy(double u, double v){
		return value(u)*value_dd(v);
	}
	static double value_dxdx_dydy(double u, double v){
		return value_dd(u)*value(v) + value(u)*value_dd(v);
	}
};


class DiscreteBspline3 : public DiscreteFilter{
public:
	DiscreteBspline3(){
		radius_h = 1;
		radius_w = 1;
	}
	double value(int r){
		if (r == -1){
			 return 1.f/6.f;
		}
		if (r == 0){
			return 2.f / 3.f;
		}
		if (r == 1){
			return 1.f / 6.f;
		}
		return 0.f;
	}
	double value(int u, int v){
		return value(u)*value(v);
	}
};

#endif //FILTER_INCLUDED