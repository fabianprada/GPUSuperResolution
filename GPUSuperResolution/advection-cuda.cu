/*
Copyright (c) 2018, Fabian Prada
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution.

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/


#include "advection-cuda.cuh"

////////////////////////////////////////////////// FILTERING //////////////////////////////////////////////////

#define BSPLINE5 0
#define BSPLINE3 0
#define GAUSSIAN 1

#if GAUSSIAN
#define PI_2 6.28318530718
#define PI_2_SQRT 2.50662827463

#define SIGMA_SQUARED 2
#define SIGMA  2
#define FILTER_RADIUS_I 6
#define FILTER_RADIUS_F 6
__constant__ float GaussianInteger[2*FILTER_RADIUS_I];

#elif BSPLINE3
#define FILTER_RADIUS_I 2
#define FILTER_RADIUS_F 2.0
#elif BSPLINE5
#define FILTER_RADIUS_I 3
#define FILTER_RADIUS_F 3.0
#else
#define FILTER_RADIUS_I 1
#define FILTER_RADIUS_F 1.0
#endif

__device__ float BilinearFilter(const float x)
{
	float r = 1.f - abs(x);
	if(r > 0.f)
		return r;
	else
		return 0.f;
}

__device__ float BilinearFilter(const float x, const float y)
{
   return BilinearFilter(x)*BilinearFilter(y);
}

__device__ float Bspline3(const float x)
{
	float r = abs(x);
	if (r < 1.f) return (4.f + r*r*(-6.f + 3.f*r)) / 6.f;
	else if (r < 2.f) return  (8.f + r*(-12.f + (6.f - r)*r)) / 6.f;
	else return 0.f;
}

__device__ float Bspline3_d(const float x)
{
		float sign_x = x > 0.f ? 1.f : -1.f;
		float r = abs(x);
		if (r < 1.f){
			return (r*(r*3.f - 4.f) / 2.f)*sign_x;
		}
		else if (r < 2.f){
			return ((r*(4.f - r) - 4.f) / 2.f)*sign_x;
		}
		else return 0.f;
}

__device__ float Bspline3_dd(const float x)
{
	float r = abs(x);
	if (r < 1.f){
		return (3.f*r-2.f);
	}
	else if (r < 2.f){
		return 2.f -r;
	}
	else return 0.f;
}

__device__ float Bspline3(const float x, const float y)
{
	return Bspline3(x)*Bspline3(y);
}


__device__ float Bspline5(const float x)
{
	float r = abs(x);
	if (r < 1.f) return (66.f + r*r*(-60.f
		+ (30.f - 10.f*r)*r*r)) / 120.f;
	else if (r < 2.f) return (51.f + r*(75.f + r*(-210.f
		+ r*(150.f + r*(-45.f + 5.f*r))))) / 120.f;
	else if (r < 3.f) return (243.f + r*(-405.f + r*(270.f
		+ r*(-90.f + (15.f - r)*r)))) / 120.f;
	else return 0.f;
}

__device__ float Bspline5_d(const float x)
{
	float sign_r = x > 0.f ? 1.f : -1.f;
	float r = abs(x);
	if (r < 1.f){
		return  (r*(-1.f + (r*r*(1.f - (5.f*r/ 12.f)))))*sign_r;
	}
	else if (r < 2.f){
		return  ((15.f + r*(-84.f + r*(90.f + r*(-36.f + 5.f*r))))/24.f)*sign_r;
	}
	else if (r < 3.f) return (-(r - 3.f)*(r - 3.f)*(r - 3.f)*(r - 3.f)/ 24.f)*sign_r;
	else return 0.f;
}

__device__ float Bspline5_dd(const float x)
{
	float r = abs(x);
	if (r < 1.f){
		return  -1.f + r*r*(3.f  - (5.f*r/ 3.f));
	}
	else if (r < 2.f){
		return  (-21.f + r*(45 + r*(-27.f + 5*r)))/6.f;
	}
	else if (r < 3.f) return  (-(r - 3.f)*(r - 3.f)*(r - 3.f) / 6.f);
	else return 0.f;
}

#if GAUSSIAN
__device__ float Gaussian(const float u ,const float v)
{
    return exp(-(u*u + v*v) / (2.f*SIGMA_SQUARED)) / (PI_2*SIGMA_SQUARED);
}
__device__ float Gaussian(const float u)
{
    return exp(-(u*u) / (2.f*SIGMA_SQUARED)) / (PI_2_SQRT*SIGMA);
}

extern "C"
void SetGaussianInteger()
{
    float values[2*FILTER_RADIUS_I];

    for (int i = -FILTER_RADIUS_I + 1; i < FILTER_RADIUS_I + 1; i++){
       values[i+FILTER_RADIUS_I-1] = exp(-((float)(i*i)) / (2.f*SIGMA_SQUARED));
    }
	cutilSafeCall(cudaMemcpyToSymbol(GaussianInteger, values, sizeof(float)*(2*FILTER_RADIUS_I)));
}

#endif


////////////////////////////////////////////////// CASTING //////////////////////////////////////////////////

__device__  void AdditionFloat4 (float4 & in,const float4 a,const float4 b)
{
	in.x = a.x +b.x;
	in.y = a.y +b.y;
	in.z = a.z +b.z;
	in.w = a.w +b.w;
}

__device__  float3 SubstractionFloat4ToFloat3(const float4 a,const float4 b)
{
	return make_float3(a.x - b.x,a.y - b.y, a.z - b.z);
}

__device__  float4 SubstractionFloat4(const float4 a,const float4 b)
{
	return make_float4(a.x - b.x,a.y - b.y, a.z - b.z,a.w - b.w);
}

__device__  float3 ScalarMultiplicationFloat3(const float3 a,const float s)
{
	return make_float3(a.x*s,a.y*s,a.z*s);
}

__device__ void NormalFloat2UChar(const float x, unsigned char & u)
{
	int i = (int)floor(x*255.f);
	if (i < 0)
		i = 0;
	else if (i >255)
		i = 255;

	u = (unsigned char) i;
}

__device__ unsigned char NormalFloat2UChar(const float x)
{
	int i = (int)floor(x*255.f);
	if (i < 0)
		i = 0;
	else if (i >255)
		i = 255;

	return (unsigned char) i;
}



////////////////////////////////////////////////// COMMOM METHODS //////////////////////////////////////////////////

__global__ void SamplePreviousLevel_Kernel(float4 * color_field_buffer,float * scalar_field_buffer,const float corner_w,const float corner_h, const float relative_scale, const int p_width, const int p_height)
{
	int global_id_x= blockIdx.x*BLOCK_SIZE + threadIdx.x;
	int global_id_y= blockIdx.y*BLOCK_SIZE + threadIdx.y;

	if( global_id_x < p_width && global_id_y < p_height ){

		float pos_w = corner_w + relative_scale*((float)(global_id_x)+0.5f);
		float pos_h = corner_h + relative_scale*((float)(global_id_y)+0.5f);

		float4 sample =  tex2D(color_texture,pos_w,pos_h);
		int write_pos = global_id_x + global_id_y*p_width;

		color_field_buffer[write_pos] = sample;
		scalar_field_buffer[write_pos] =sample.x*0.3f + sample.y*0.59f + sample.z*0.11f; 
	}
}


__global__ void SampleGlobalAdvectionField_Kernel( float2* dst, float * norms, const unsigned int imgWidth,const unsigned int imgHeight, const float finv_imgWidht,  const float finv_imgHeight, const float corner_w,const float corner_h,const float scale,  bool p_normalize_advection_field )
{
	if( (blockIdx.x*BLOCK_SIZE + threadIdx.x) < imgWidth && (blockIdx.y*BLOCK_SIZE + threadIdx.y) < imgHeight ){

	float local_pos_w = (float)(blockIdx.x*BLOCK_SIZE + threadIdx.x) + 0.5f;
	float local_pos_h = (float)(blockIdx.y*BLOCK_SIZE + threadIdx.y) + 0.5f;

	float global_pos_w = corner_w + scale*local_pos_w;
	float global_pos_h = corner_h + scale*local_pos_h;
	
	float cw = floor(global_pos_w);
	float ch = floor(global_pos_h);

	float dw = global_pos_w-cw;
	float dh = global_pos_h-ch;

	float filter_pos_w;
	float filter_pos_h;

	float grad_x = 0.f; 
	float grad_y = 0.f;
	float hessian_xx = 0.f;
	float hessian_yy = 0.f;
	float hessian_xy = 0.f;

	float filter_value;
	float scalar_value;

#if BSPLINE3
	float bspline_x;
	float bspline_y;
	float bspline_d_x;
	float bspline_d_y;
#endif

#if BSPLINE5
	float bspline_x;
	float bspline_y;
	float bspline_d_x;
	float bspline_d_y;
#endif

	for (int ih = -FILTER_RADIUS_I + 1; ih < FILTER_RADIUS_I + 1; ih++){
			for (int iw = -FILTER_RADIUS_I + 1; iw < FILTER_RADIUS_I + 1; iw++){
				filter_pos_w = dw - (float)(iw); // CONVOLUTION POSITION
				filter_pos_h = dh - (float)(ih);
				scalar_value = tex2D(scalar_texture,cw + (float)iw,ch + (float)ih);
#if GAUSSIAN
				filter_value =  Gaussian(filter_pos_w,filter_pos_h);
				grad_x +=(-filter_pos_w / SIGMA_SQUARED)*filter_value*scalar_value;
				grad_y +=(-filter_pos_h / SIGMA_SQUARED)*filter_value*scalar_value;
				hessian_xx += ((filter_pos_w*filter_pos_w - SIGMA_SQUARED)/(SIGMA_SQUARED*SIGMA_SQUARED))*filter_value*scalar_value;
				hessian_yy += ((filter_pos_h*filter_pos_h - SIGMA_SQUARED)/(SIGMA_SQUARED*SIGMA_SQUARED))*filter_value*scalar_value;
				hessian_xy += ((filter_pos_w*filter_pos_h)/(SIGMA_SQUARED*SIGMA_SQUARED))*filter_value*scalar_value;
#elif BSPLINE3
				bspline_x =Bspline3(filter_pos_w);
				bspline_y =Bspline3(filter_pos_h);
				bspline_d_x = Bspline3_d(filter_pos_w);
				bspline_d_y = Bspline3_d(filter_pos_h);
				grad_x += bspline_d_x*bspline_y*scalar_value;
				grad_y += bspline_d_y*bspline_x*scalar_value;
				hessian_xx += Bspline3_dd(filter_pos_w)*bspline_y*scalar_value;
				hessian_yy += Bspline3_dd(filter_pos_h)*bspline_x*scalar_value;
				hessian_xy += bspline_d_x*bspline_d_y*scalar_value;
#elif BSPLINE5
				bspline_x =Bspline5(filter_pos_w);
				bspline_y =Bspline5(filter_pos_h);
				bspline_d_x = Bspline5_d(filter_pos_w);
				bspline_d_y = Bspline5_d(filter_pos_h);
				grad_x += bspline_d_x*bspline_y*scalar_value;
				grad_y += bspline_d_y*bspline_x*scalar_value;
				hessian_xx += Bspline5_dd(filter_pos_w)*bspline_y*scalar_value;
				hessian_yy += Bspline5_dd(filter_pos_h)*bspline_x*scalar_value;
				hessian_xy += bspline_d_x*bspline_d_y*scalar_value;
#endif
			}
	}
	scalar_value =hessian_xx*grad_x*grad_x + 2.f*hessian_xy*grad_x*grad_y +hessian_yy*grad_y*grad_y;
	grad_x *= scalar_value;
	grad_y *= scalar_value;

	float grad_norm = sqrt(grad_x*grad_x + grad_y*grad_y);

	if(p_normalize_advection_field && grad_norm > 0.f){
			grad_x/=grad_norm;
			grad_y/=grad_norm;
			norms[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth] = 1.f;
	}
	else{
		norms[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth] = grad_norm;
	}

	dst[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth].x  = grad_x;
	dst[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth].y  = grad_y;
	}
}

__global__ void SampleGlobalAdvectionField_Constant_Kernel( float2* dst, float * norms, const unsigned int imgWidth,const unsigned int imgHeight, const float finv_imgWidht,  const float finv_imgHeight, const float corner_w,const float corner_h,const float scale,  bool p_normalize_advection_field )
{
	if( (blockIdx.x*BLOCK_SIZE + threadIdx.x) < imgWidth && (blockIdx.y*BLOCK_SIZE + threadIdx.y) < imgHeight ){

	float dw = corner_w + scale*((float)(blockIdx.x*BLOCK_SIZE + threadIdx.x) + 0.5f);
	float dh = corner_h + scale*((float)(blockIdx.y*BLOCK_SIZE + threadIdx.y) + 0.5f);
	
	float cw = floor(dw);
	float ch = floor(dh);

	dw -= cw;
	dh -= ch;

	float filter_pos_w;
	float filter_pos_h;

	float exp_2dw = exp(2.f*dw / (2.f*SIGMA_SQUARED));
	float exp_2dh = exp(2.f*dh / (2.f*SIGMA_SQUARED));

	float exp_dw_first = exp(-(dw*(dw+2.f*(float)(FILTER_RADIUS_I-1)))/(2.f*SIGMA_SQUARED))/(PI_2*SIGMA_SQUARED);
	float exp_dh_first = exp(-(dh*(dh+2.f*(float)(FILTER_RADIUS_I-1)))/(2.f*SIGMA_SQUARED));

	float grad_x = 0.f; 
	float grad_y = 0.f;
	float hessian_xx = 0.f;
	float hessian_yy = 0.f;
	float hessian_xy = 0.f;

	float filter_value;
	float temporal_filter_value;
	float scalar_value;

	for (int ih = -FILTER_RADIUS_I + 1; ih < FILTER_RADIUS_I + 1; ih++){
		    temporal_filter_value = exp_dw_first*exp_dh_first*GaussianInteger[ih + FILTER_RADIUS_I - 1];
			filter_pos_h = dh - (float)(ih);
			for (int iw = -FILTER_RADIUS_I + 1; iw < FILTER_RADIUS_I + 1; iw++){

				filter_pos_w = dw - (float)(iw); // CONVOLUTION POSITION
				scalar_value = tex2D(scalar_texture,cw + (float)iw,ch + (float)ih);

				filter_value = temporal_filter_value * GaussianInteger[iw + FILTER_RADIUS_I - 1];
				grad_x +=(-filter_pos_w / SIGMA_SQUARED)*filter_value*scalar_value;
				grad_y +=(-filter_pos_h / SIGMA_SQUARED)*filter_value*scalar_value;
				hessian_xx += ((filter_pos_w*filter_pos_w - SIGMA_SQUARED)/(SIGMA_SQUARED*SIGMA_SQUARED))*filter_value*scalar_value;
				hessian_yy += ((filter_pos_h*filter_pos_h - SIGMA_SQUARED)/(SIGMA_SQUARED*SIGMA_SQUARED))*filter_value*scalar_value;
				hessian_xy += ((filter_pos_w*filter_pos_h)/(SIGMA_SQUARED*SIGMA_SQUARED))*filter_value*scalar_value;

				temporal_filter_value*=exp_2dw;
			}
		    exp_dh_first*=exp_2dh;
	}

	scalar_value =hessian_xx*grad_x*grad_x + 2.f*hessian_xy*grad_x*grad_y +hessian_yy*grad_y*grad_y;
	grad_x *= scalar_value;
	grad_y *= scalar_value;

	float grad_norm = sqrt(grad_x*grad_x + grad_y*grad_y);

	if(p_normalize_advection_field && grad_norm > 0.f){
			grad_x/=grad_norm;
			grad_y/=grad_norm;
			norms[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth] = 1.f;
	}
	else{
		norms[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth] = grad_norm;
	}

	dst[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth].x  = grad_x;
	dst[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth].y  = grad_y;
	}
}


__global__ void SampleGlobalMeanGradientField_Kernel( float2* dst, float * norms, const unsigned int imgWidth,const unsigned int imgHeight, const float finv_imgWidht,  const float finv_imgHeight, const float corner_w,const float corner_h,const float scale,  bool p_normalize_advection_field )
{
	if( (blockIdx.x*BLOCK_SIZE + threadIdx.x) < imgWidth && (blockIdx.y*BLOCK_SIZE + threadIdx.y) < imgHeight ){

	float local_pos_w = (float)(blockIdx.x*BLOCK_SIZE + threadIdx.x) + 0.5f;
	float local_pos_h = (float)(blockIdx.y*BLOCK_SIZE + threadIdx.y) + 0.5f;

	float global_pos_w = corner_w + scale*local_pos_w;
	float global_pos_h = corner_h + scale*local_pos_h;
	
	float cw = floor(global_pos_w);
	float ch = floor(global_pos_h);

	float dw = global_pos_w-cw;
	float dh = global_pos_h-ch;

	float filter_pos_w;
	float filter_pos_h;

	float grad_x = 0.f; 
	float grad_y = 0.f;

	float filter_value;
	float scalar_value;


	for (int iw = -FILTER_RADIUS_I + 1; iw < FILTER_RADIUS_I + 1; iw++){
			for (int ih = -FILTER_RADIUS_I + 1; ih < FILTER_RADIUS_I + 1; ih++){
				filter_pos_w = dw - (float)(iw); // CONVOLUTION POSITION
				filter_pos_h = dh - (float)(ih);
				scalar_value = tex2D(scalar_texture,cw + (float)iw,ch + (float)ih);
#if GAUSSIAN
				filter_value =  Gaussian(filter_pos_w,filter_pos_h);
				grad_x +=(-filter_pos_w / SIGMA_SQUARED)*filter_value*scalar_value;
				grad_y +=(-filter_pos_h / SIGMA_SQUARED)*filter_value*scalar_value;
#elif BSPLINE3
				grad_x += Bspline3_d(filter_pos_w)*Bspline3(filter_pos_h)*scalar_value;
				grad_y += Bspline3_d(filter_pos_h)*Bspline3(filter_pos_w)*scalar_value;
#elif BSPLINE5
				grad_x += Bspline5_d(filter_pos_w)*Bspline5(filter_pos_h)*scalar_value;
				grad_y += Bspline5_d(filter_pos_h)*Bspline5(filter_pos_w)*scalar_value;
#else
#endif
			}
	}


	float grad_norm = sqrt(grad_x*grad_x + grad_y*grad_y);

	if(p_normalize_advection_field && grad_norm > 0.f){
			grad_x/=grad_norm;
			grad_y/=grad_norm;
			norms[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth] = 1.f;
	}
	else{
		norms[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth] = grad_norm;
	}

	dst[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth].x  = grad_x;
	dst[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth].y  = grad_y;
	}
}

__global__ void SampleLocalAdvectionField_Kernel( float2* dst, float * norms, const unsigned int imgWidth,const unsigned int imgHeight, bool p_normalize_advection_field)
{
	if( (blockIdx.x*BLOCK_SIZE + threadIdx.x) < imgWidth && (blockIdx.y*BLOCK_SIZE + threadIdx.y) < imgHeight ){

	float local_pos_w = (float)(blockIdx.x*BLOCK_SIZE + threadIdx.x) + 0.5f;
	float local_pos_h = (float)(blockIdx.y*BLOCK_SIZE + threadIdx.y) + 0.5f;
	
	float filter_pos_w;
	float filter_pos_h;

	float grad_x = 0.f; 
	float grad_y = 0.f;
	float hessian_xx = 0.f;
	float hessian_yy = 0.f;
	float hessian_xy = 0.f;

	float filter_value;
	float scalar_value;

	for (int ih = -FILTER_RADIUS_I + 1; ih < FILTER_RADIUS_I + 1; ih++){
			for (int iw = -FILTER_RADIUS_I + 1; iw < FILTER_RADIUS_I + 1; iw++){
				filter_pos_w = -(float)(iw); // CONVOLUTION POSITION
				filter_pos_h = -(float)(ih);
				scalar_value = tex2D(scalar_texture,local_pos_w - filter_pos_w,local_pos_h - filter_pos_h);
#if GAUSSIAN
				filter_value =  Gaussian(filter_pos_w,filter_pos_h);
				grad_x -=(filter_pos_w / SIGMA_SQUARED)*filter_value*scalar_value;
				grad_y -=(filter_pos_h / SIGMA_SQUARED)*filter_value*scalar_value;
				hessian_xx += ((filter_pos_w*filter_pos_w - SIGMA_SQUARED)/(SIGMA_SQUARED*SIGMA_SQUARED))*filter_value*scalar_value;
				hessian_yy += ((filter_pos_h*filter_pos_h - SIGMA_SQUARED)/(SIGMA_SQUARED*SIGMA_SQUARED))*filter_value*scalar_value;
				hessian_xy += ((filter_pos_w*filter_pos_h)/(SIGMA_SQUARED*SIGMA_SQUARED))*filter_value*scalar_value;
#elif BSPLINE3
#else
#endif
			}
	}
	scalar_value =hessian_xx*grad_x*grad_x + 2.f*hessian_xy*grad_x*grad_y +hessian_yy*grad_y*grad_y;
	grad_x *= scalar_value;
	grad_y *= scalar_value;

	float grad_norm = sqrt(grad_x*grad_x + grad_y*grad_y);

	if(p_normalize_advection_field && grad_norm > 0.f){
			grad_x/=grad_norm;
			grad_y/=grad_norm;
			norms[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth] = 1.f;
	}
	else{
		norms[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth] = grad_norm;
	}

	dst[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth].x  = grad_x;
	dst[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth].y  = grad_y;
	}
}

__global__ void SampleLocalMeanGradientField_Kernel( float2* dst, float * norms, const unsigned int imgWidth,const unsigned int imgHeight, bool p_normalize_advection_field)
{
	if( (blockIdx.x*BLOCK_SIZE + threadIdx.x) < imgWidth && (blockIdx.y*BLOCK_SIZE + threadIdx.y) < imgHeight ){

	float local_pos_w = (float)(blockIdx.x*BLOCK_SIZE + threadIdx.x) + 0.5f;
	float local_pos_h = (float)(blockIdx.y*BLOCK_SIZE + threadIdx.y) + 0.5f;
	
	float filter_pos_w;
	float filter_pos_h;

	float grad_x = 0.f; 
	float grad_y = 0.f;
	float filter_value;
	float scalar_value;

	for (int ih = -FILTER_RADIUS_I + 1; ih < FILTER_RADIUS_I + 1; ih++){
			for (int iw = -FILTER_RADIUS_I + 1; iw < FILTER_RADIUS_I + 1; iw++){
				filter_pos_w = -(float)(iw); // CONVOLUTION POSITION
				filter_pos_h = -(float)(ih);
				scalar_value = tex2D(scalar_texture,local_pos_w - filter_pos_w,local_pos_h -filter_pos_h);
#if GAUSSIAN
				filter_value =  Gaussian(filter_pos_w,filter_pos_h);
				grad_x -=(filter_pos_w/ SIGMA_SQUARED)*filter_value*scalar_value;
				grad_y -=(filter_pos_h / SIGMA_SQUARED)*filter_value*scalar_value;
#elif BSPLINE3
				grad_x +=Bspline3_d(filter_pos_w)*Bspline3(filter_pos_h)*scalar_value;
				grad_y +=Bspline3_d(filter_pos_h)*Bspline3(filter_pos_w)*scalar_value;
#elif BSPLINE5
				grad_x += Bspline5_d(filter_pos_w)*Bspline5(filter_pos_h)*scalar_value;
				grad_y += Bspline5_d(filter_pos_h)*Bspline5(filter_pos_w)*scalar_value;
#else
#endif
			}
	}

	float grad_norm = sqrt(grad_x*grad_x + grad_y*grad_y);

	if(p_normalize_advection_field && grad_norm > 0.f){
			grad_x/=grad_norm;
			grad_y/=grad_norm;
			norms[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth] = 1.f;
	}
	else{
		norms[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth] = grad_norm;
	}

	dst[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth].x  = grad_x;
	dst[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth].y  = grad_y;
	}
}

__global__ void SampleLocalMeanGradientFieldLinear_Kernel( float2* dst, float * norms, const unsigned int imgWidth,const unsigned int imgHeight, bool p_normalize_advection_field)
{
	if( (blockIdx.x*BLOCK_SIZE + threadIdx.x) < imgWidth && (blockIdx.y*BLOCK_SIZE + threadIdx.y) < imgHeight ){

	float local_pos_w = (float)(blockIdx.x*BLOCK_SIZE + threadIdx.x) + 0.5f;
	float local_pos_h = (float)(blockIdx.y*BLOCK_SIZE + threadIdx.y) + 0.5f;
	
	float grad_x = tex2D(scalar_texture,local_pos_w + 1.f,local_pos_h) - tex2D(scalar_texture,local_pos_w - 1.f,local_pos_h); 
	float grad_y = tex2D(scalar_texture,local_pos_w ,local_pos_h + 1.f) - tex2D(scalar_texture,local_pos_w,local_pos_h - 1.f); 

	float grad_norm = sqrt(grad_x*grad_x + grad_y*grad_y);

	if(p_normalize_advection_field && grad_norm > 0.f){
			grad_x/=grad_norm;
			grad_y/=grad_norm;
			norms[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth] = 1.f;
	}
	else{
		norms[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth] = grad_norm;
	}

	dst[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth].x  = grad_x;
	dst[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth].y  = grad_y;
	}
}

__global__ void SetTexture_Kernel( unsigned char * uchar4_buffer_d, float * buffer_r_d,float * buffer_g_d,float * buffer_b_d, const unsigned int imgWidth,const unsigned int imgHeight)
{
	int global_id_x= blockIdx.x*BLOCK_SIZE + threadIdx.x;
	int global_id_y= blockIdx.y*BLOCK_SIZE + threadIdx.y;

	if( global_id_x < imgWidth && global_id_y < imgHeight ){

	unsigned int index = global_id_x + global_id_y*imgWidth;
	unsigned int writting_index = 4*index;
	uchar4_buffer_d[writting_index] = NormalFloat2UChar(buffer_r_d[index]);
	uchar4_buffer_d[writting_index + 1] = NormalFloat2UChar(buffer_g_d[index]);
	uchar4_buffer_d[writting_index + 2] = NormalFloat2UChar(buffer_b_d[index]);
	uchar4_buffer_d[writting_index + 3] = 255;
	}

	//if( (blockIdx.x*BLOCK_SIZE + threadIdx.x) < imgWidth && (blockIdx.y*BLOCK_SIZE + threadIdx.y) < imgHeight ){
	//	float pos_w = corner_w + scale*((float)(blockIdx.x*BLOCK_SIZE + threadIdx.x) + 0.5f)*finv_imgWidht;
	//	float pos_h = corner_h + scale*((float)(blockIdx.y*BLOCK_SIZE + threadIdx.y) + 0.5f)*finv_imgHeight;
	//	dst[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth] = tex2D(scalar_texture,pos_w,pos_h);
	//}
}


__global__ void FloatToUchar4( float * src, unsigned char * dst,const unsigned int imgWidth, const unsigned int imgHeight)
{
	if( (blockIdx.x*BLOCK_SIZE + threadIdx.x) < imgWidth && (blockIdx.y*BLOCK_SIZE + threadIdx.y) < imgHeight ){
	
	int position =(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth; 
	int i_color = (int)floor(src[position]*255.f);
	if (i_color < 0)
		i_color = 0;
	else if (i_color>255)
		i_color = 255;

    unsigned char value = (unsigned char)(i_color);
	dst[4*position] = value;
	dst[4*position + 1] = value;
	dst[4*position + 2] = value;
	dst[4*position + 3] = 255;
	}
}

__global__ void BackwardSamplingFromGlobal_Kernel(float4 * dst, const unsigned int imgWidth, const unsigned int imgHeight, const float corner_w,const float corner_h,const float scale ,const unsigned int num_iter, const unsigned int supersampling, const float angle_tolerance , const float step_amplification)
{
	//int global_id_x= blockIdx.x*BLOCK_SIZE + threadIdx.x;
	//int global_id_y= blockIdx.y*BLOCK_SIZE + threadIdx.y;

	if( (blockIdx.x*BLOCK_SIZE + threadIdx.x) < imgWidth && (blockIdx.y*BLOCK_SIZE + threadIdx.y) < imgHeight ){

	float local_pos_w = (float)(blockIdx.x*BLOCK_SIZE + threadIdx.x) + 0.5f;
	float local_pos_h = (float)(blockIdx.y*BLOCK_SIZE + threadIdx.y) + 0.5f;

	float4 cummulative_color = make_float4(0.f,0.f,0.f,0.f);
	float4 sampled_color;

	float cummulative_weight = 0.f;
	float weight;

	float pos_w;
	float pos_h;

	float dw;
	float dh;

	float2 advection_direction;
	float2 last_advection_direction;

	float angle2;
	float advection_direction_norm2;
	float last_advection_direction_norm2;

	for( int iw =0; iw<supersampling; iw++){
		for( int jw =0; jw<supersampling; jw++){

	//dw = - FILTER_RADIUS_F + 2.f*FILTER_RADIUS_F*((float)iw + 0.5f)/((float)supersampling);
	//dh = - FILTER_RADIUS_F + 2.f*FILTER_RADIUS_F*((float)jw + 0.5f)/((float)supersampling);

    dw = -1.f + 2.f*((float)iw + 0.5f)/((float)supersampling);
	dh = -1.f + 2.f*((float)jw + 0.5f)/((float)supersampling);

	pos_w = local_pos_w + dw;
	pos_h = local_pos_h + dh;
	
	advection_direction = tex2D(advection_texture,pos_w,pos_h);
	advection_direction_norm2 = advection_direction.x*advection_direction.x + advection_direction.y*advection_direction.y;

	angle2 = 1.f;

	for(int i =0; i<num_iter; i++){
		if(angle2>angle_tolerance){
			float potential_pos_w = pos_w - advection_direction.x*step_amplification;
			float potential_pos_h = pos_h - advection_direction.y*step_amplification;
			
			if((potential_pos_w > 0.f && potential_pos_w < (float)(imgWidth-1)) && (potential_pos_h > 0.f && potential_pos_h < (float)(imgHeight-1))){
				
			pos_w = potential_pos_w;
			pos_h = potential_pos_h;

			last_advection_direction = advection_direction;
			last_advection_direction_norm2 = advection_direction_norm2;

			advection_direction = tex2D(advection_texture,pos_w,pos_h);
			advection_direction_norm2 = advection_direction.x*advection_direction.x + advection_direction.y*advection_direction.y;

			if(advection_direction_norm2>0.f){
					angle2 = last_advection_direction.x*advection_direction.x + last_advection_direction.y*advection_direction.y;
					angle2 = angle2>0.f ? angle2*angle2 : -angle2*angle2;
					angle2 /=(last_advection_direction_norm2*advection_direction_norm2);
				}
				else{
					angle2 = angle_tolerance - 1.f;
				}
		    }
	   }
	}

	//pos_w = corner_w + scale*(pos_w+0.5f)*finv_imgWidht;
	//pos_h = corner_h + scale*(pos_h+0.5f)*finv_imgHeight;
	pos_w = corner_w + scale*(pos_w); // Unnormalized color texture
	pos_h = corner_h + scale*(pos_h);
	sampled_color = tex2D(color_texture,pos_w,pos_h);

//#if GAUSSIAN
//	weight = Gaussian(dw,dh);
//#elif BSPLINE3
//	weight = Bspline3(dw,dh);
//#else
//	weight = BilinearFilter(dw,dh);
//#endif

	weight = BilinearFilter(dw,dh);

	cummulative_weight += weight; 
	cummulative_color.x += sampled_color.x*weight;
	cummulative_color.y += sampled_color.y*weight;
	cummulative_color.z += sampled_color.z*weight;

		}
	}

	cummulative_color.x/=cummulative_weight;
	cummulative_color.y/=cummulative_weight;
	cummulative_color.z/=cummulative_weight;

	int position =(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth; 

	dst[position] = cummulative_color;
	}
}

__global__ void BackwardSamplingFromLocal_Kernel(float4 * dst, const unsigned int imgWidth, const unsigned int imgHeight,const unsigned int num_iter, const unsigned int supersampling , const float angle_tolerance,const float step_amplification)
{

	if( (blockIdx.x*BLOCK_SIZE + threadIdx.x) < imgWidth && (blockIdx.y*BLOCK_SIZE + threadIdx.y) < imgHeight ){

	float local_pos_w = (float)(blockIdx.x*BLOCK_SIZE + threadIdx.x) + 0.5f;
	float local_pos_h = (float)(blockIdx.y*BLOCK_SIZE + threadIdx.y) + 0.5f;

	float4 cummulative_color = make_float4(0.f,0.f,0.f,0.f);
	float4 sampled_color;

	float cummulative_weight = 0.f;
	float weight;

	float pos_w;
	float pos_h;

	float dw;
	float dh;

	float2 advection_direction;
	float2 last_advection_direction;

	float angle2;
	float advection_direction_norm2;
	float last_advection_direction_norm2;

	for( int iw =0; iw<supersampling; iw++){
		for( int jw =0; jw<supersampling; jw++){

	//dw = - FILTER_RADIUS_F + 2.f*FILTER_RADIUS_F*((float)iw + 0.5f)/((float)supersampling);
	//dh = - FILTER_RADIUS_F + 2.f*FILTER_RADIUS_F*((float)jw + 0.5f)/((float)supersampling);

    dw = -1.f + 2.f*((float)iw + 0.5f)/((float)supersampling);
	dh = -1.f + 2.f*((float)jw + 0.5f)/((float)supersampling);

	pos_w = local_pos_w + dw;
	pos_h = local_pos_h + dh;
	
	advection_direction = tex2D(advection_texture,pos_w,pos_h);
	advection_direction_norm2 = advection_direction.x*advection_direction.x + advection_direction.y*advection_direction.y;

	angle2 = 1.f;

	for(int i=0; i<num_iter; i++){		
		if(angle2>angle_tolerance){
			float potential_pos_w = pos_w - advection_direction.x*step_amplification;
			float potential_pos_h = pos_h - advection_direction.y*step_amplification;
			if((potential_pos_w > 0.f && potential_pos_w < (float)(imgWidth-1)) && (potential_pos_h > 0.f && potential_pos_h < (float)(imgHeight-1))){
				
				pos_w = potential_pos_w;
				pos_h = potential_pos_h;

				last_advection_direction = advection_direction;
				last_advection_direction_norm2 = advection_direction_norm2;

				advection_direction = tex2D(advection_texture,pos_w,pos_h);
				advection_direction_norm2 = advection_direction.x*advection_direction.x + advection_direction.y*advection_direction.y;
				if(advection_direction_norm2>0.f){
					angle2 = last_advection_direction.x*advection_direction.x + last_advection_direction.y*advection_direction.y;
					angle2 = angle2>0.f ? angle2*angle2 : -angle2*angle2;
					angle2 /=(last_advection_direction_norm2*advection_direction_norm2);
				}
				else{
					angle2 = angle_tolerance - 1.f;
				}
			}
		}
	}

	sampled_color = tex2D(color_texture,pos_w,pos_h);

//#if GAUSSIAN
//	weight = Gaussian(dw,dh);
//#elif BSPLINE3
//	weight = Bspline3(dw,dh);
//#else
//	weight = BilinearFilter(dw,dh);
//#endif
	weight = BilinearFilter(dw,dh);

	cummulative_weight += weight; 
	cummulative_color.x += sampled_color.x*weight;
	cummulative_color.y += sampled_color.y*weight;
	cummulative_color.z += sampled_color.z*weight;

		}
	}

	cummulative_color.x/=cummulative_weight;
	cummulative_color.y/=cummulative_weight;
	cummulative_color.z/=cummulative_weight;

	int position =(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth; 

	dst[position] = cummulative_color;
	}
}

///////////////////////////////////////// VISUALIZATION TOOLS ///////////////////////////////////////// 

#include "gltools.h"

__global__ void SetVectorField_Kernel( float * dst, float2 * vector_field, const unsigned int imgWidth,const unsigned int imgHeight, const float finv_imgWidht, const float finv_imgHeight, bool normalize_for_visualization, bool clamp_for_visualization, const float clamping_threshold, const float visual_amplification)
{
	int global_id_x= blockIdx.x*BLOCK_SIZE + threadIdx.x;
	int global_id_y= blockIdx.y*BLOCK_SIZE + threadIdx.y;

	if( global_id_x < imgWidth && global_id_y < imgHeight ){

	float unnormalized_pos_w = ((float)(global_id_x) + 0.5f);
	float unnormalized_pos_h = ((float)(global_id_y) + 0.5f);

	unsigned int reading_index = global_id_x + global_id_y*imgWidth;
	unsigned int writting_index = 9*reading_index;
	
	float2 vector_direction = vector_field[reading_index];


	float vector_direction_norm = vector_direction.x*vector_direction.x + vector_direction.y*vector_direction.y;
	if(vector_direction_norm > 0.f){
		vector_direction_norm = sqrt(vector_direction_norm);
		if(normalize_for_visualization){
			vector_direction.x*= (visual_amplification/vector_direction_norm);
			vector_direction.y*= (visual_amplification/vector_direction_norm);
		}
		else if(clamp_for_visualization && vector_direction_norm > clamping_threshold){
			vector_direction.x*= (visual_amplification*clamping_threshold/vector_direction_norm);
			vector_direction.y*= (visual_amplification*clamping_threshold/vector_direction_norm);
		}
		else{
			vector_direction.x*= visual_amplification;
			vector_direction.y*= visual_amplification;
		}
	}

	dst[writting_index] = (unnormalized_pos_w + vector_direction.x)*finv_imgWidht;
	dst[writting_index + 1] = (unnormalized_pos_h + vector_direction.y)*finv_imgHeight;
	dst[writting_index + 2] = 0.f;
	dst[writting_index + 3] = (unnormalized_pos_w - vector_direction.y * 0.25f)*finv_imgWidht;
	dst[writting_index + 4] = (unnormalized_pos_h + vector_direction.x * 0.25f)*finv_imgHeight;
	dst[writting_index + 5] = 0.f;
	dst[writting_index + 6] = (unnormalized_pos_w + vector_direction.y * 0.25f)*finv_imgWidht;
	dst[writting_index + 7] = (unnormalized_pos_h - vector_direction.x * 0.25f)*finv_imgHeight;
	dst[writting_index + 8] = 0.f;
	}
}

void GLTools::TransformFloat2ArrayToGLVectorBuffer_CUDA(cudaGraphicsResource_t& cuda_resource, float2 * vector_buffer, const unsigned int imgWidth,const unsigned int imgHeigh, bool normalize_vector_visualization,bool clamp_for_visualization,const float clamping_threshold,const float field_amplification)
{
	cutilSafeCall( cudaGraphicsMapResources( 1 , &cuda_resource ) );
	float * vertex_buffer_handle;
	size_t num_bytes;
	cutilSafeCall( cudaGraphicsResourceGetMappedPointer((void **)&vertex_buffer_handle,&num_bytes,cuda_resource));

	unsigned int  blocksW = (unsigned int) ceilf( (float) imgWidth / (float) BLOCK_SIZE );
    unsigned int  blocksH = (unsigned int) ceilf( (float) imgHeigh /(float) BLOCK_SIZE );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );

	SetVectorField_Kernel<<< gridDim, blockDim >>>( vertex_buffer_handle, vector_buffer, imgWidth, imgHeigh, 1.f/(float)imgWidth, 1.f/(float)imgHeigh,normalize_vector_visualization,clamp_for_visualization,clamping_threshold,field_amplification);

	cutilSafeCall( cudaGraphicsUnmapResources( 1, &cuda_resource) );

}

__global__ void SetVectorField_Kernel2( float * dst, float * vector0_field, float * vector1_field, const unsigned int imgWidth,const unsigned int imgHeight, const float finv_imgWidht, const float finv_imgHeight, bool normalize_for_visualization, bool clamp_for_visualization, const float clamping_threshold, const float visual_amplification)
{
	int global_id_x= blockIdx.x*BLOCK_SIZE + threadIdx.x;
	int global_id_y= blockIdx.y*BLOCK_SIZE + threadIdx.y;

	if( global_id_x < imgWidth && global_id_y < imgHeight ){

	float unnormalized_pos_w = ((float)(global_id_x) + 0.5f);
	float unnormalized_pos_h = ((float)(global_id_y) + 0.5f);

	unsigned int reading_index = global_id_x + global_id_y*imgWidth;
	unsigned int writting_index = 9*reading_index;
	
	float2 vector_direction;
	vector_direction.x = vector0_field[reading_index];
	vector_direction.y = vector1_field[reading_index];

	float vector_direction_norm = vector_direction.x*vector_direction.x + vector_direction.y*vector_direction.y;
	if(vector_direction_norm > 0.f){
		vector_direction_norm = sqrt(vector_direction_norm);
		if(normalize_for_visualization){
			vector_direction.x*= (visual_amplification/vector_direction_norm);
			vector_direction.y*= (visual_amplification/vector_direction_norm);
		}
		else if(clamp_for_visualization && vector_direction_norm > clamping_threshold){
			vector_direction.x*= (visual_amplification*clamping_threshold/vector_direction_norm);
			vector_direction.y*= (visual_amplification*clamping_threshold/vector_direction_norm);
		}
		else{
			vector_direction.x*= visual_amplification;
			vector_direction.y*= visual_amplification;
		}
	}

	dst[writting_index] = (unnormalized_pos_w + vector_direction.x)*finv_imgWidht;
	dst[writting_index + 1] = (unnormalized_pos_h + vector_direction.y)*finv_imgHeight;
	dst[writting_index + 2] = 0.f;
	dst[writting_index + 3] = (unnormalized_pos_w - vector_direction.y * 0.25f)*finv_imgWidht;
	dst[writting_index + 4] = (unnormalized_pos_h + vector_direction.x * 0.25f)*finv_imgHeight;
	dst[writting_index + 5] = 0.f;
	dst[writting_index + 6] = (unnormalized_pos_w + vector_direction.y * 0.25f)*finv_imgWidht;
	dst[writting_index + 7] = (unnormalized_pos_h - vector_direction.x * 0.25f)*finv_imgHeight;
	dst[writting_index + 8] = 0.f;
	}
}

void GLTools::Transform2FloatArrayToGLVectorBuffer_CUDA(cudaGraphicsResource_t& cuda_resource, float * vector0_buffer,float * vector1_buffer, const unsigned int imgWidth,const unsigned int imgHeigh, bool normalize_vector_visualization,bool clamp_for_visualization,const float clamping_threshold,const float field_amplification)
{
	cutilSafeCall( cudaGraphicsMapResources( 1 , &cuda_resource ) );
	float * vertex_buffer_handle;
	size_t num_bytes;
	cutilSafeCall( cudaGraphicsResourceGetMappedPointer((void **)&vertex_buffer_handle,&num_bytes,cuda_resource));

	unsigned int  blocksW = (unsigned int) ceilf( (float) imgWidth / (float) BLOCK_SIZE );
    unsigned int  blocksH = (unsigned int) ceilf( (float) imgHeigh /(float) BLOCK_SIZE );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );

	SetVectorField_Kernel2<<< gridDim, blockDim >>>( vertex_buffer_handle, vector0_buffer,vector1_buffer,imgWidth, imgHeigh, 1.f/(float)imgWidth, 1.f/(float)imgHeigh,normalize_vector_visualization,clamp_for_visualization,clamping_threshold,field_amplification);

	cutilSafeCall( cudaGraphicsUnmapResources( 1, &cuda_resource) );
}

__global__ void SetTexture_Kernel(unsigned char * uchar4_buffer_d,float4 * float4_buffer_d,const int p_width, const int p_height){

	int global_id_x= blockIdx.x*BLOCK_SIZE + threadIdx.x;
	int global_id_y= blockIdx.y*BLOCK_SIZE + threadIdx.y;

	if( global_id_x < p_width && global_id_y < p_height ){

	unsigned int index = global_id_x + global_id_y*p_width;
	unsigned int writting_index = 4*index;
	uchar4_buffer_d[writting_index] = NormalFloat2UChar(float4_buffer_d[index].x);
	uchar4_buffer_d[writting_index + 1] = NormalFloat2UChar(float4_buffer_d[index].y);
	uchar4_buffer_d[writting_index + 2] = NormalFloat2UChar(float4_buffer_d[index].z);
	uchar4_buffer_d[writting_index + 3] = 255;
	}
}


void GLTools::TransformFloat4ArrayToGLColorBuffer_CUDA(cudaGraphicsResource_t& cuda_resource, float4 * float4_buffer, unsigned char * uchar4_buffer, const unsigned int width,const unsigned int height)
{
	cutilSafeCall( cudaGraphicsMapResources( 1, &cuda_resource ) );
	cudaArray* texture_buffer_handle;
	cutilSafeCall( cudaGraphicsSubResourceGetMappedArray( &texture_buffer_handle, cuda_resource, 0, 0 ) );

	unsigned int  blocksW = (unsigned int) ceilf( (float) width / (float) BLOCK_SIZE );
    unsigned int  blocksH = (unsigned int) ceilf( (float) height /(float) BLOCK_SIZE );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );

	SetTexture_Kernel<<< gridDim, blockDim >>>(uchar4_buffer,float4_buffer, width, height);

	cutilSafeCall( cudaMemcpyToArray( texture_buffer_handle, 0, 0, uchar4_buffer, width*height*4*sizeof(unsigned char), cudaMemcpyDeviceToDevice ) );
	cutilSafeCall( cudaGraphicsUnmapResources( 1, &cuda_resource) );
}

__global__ void CopyFloatFromFloat4_Kernel(float4 *src, float * dst, const unsigned int width,const unsigned int height, int channel)
{
	int global_id_x= blockIdx.x*BLOCK_SIZE + threadIdx.x;
	int global_id_y= blockIdx.y*BLOCK_SIZE + threadIdx.y;

	if( global_id_x < width && global_id_y < height ){
		unsigned int index = global_id_x + global_id_y*width;
		if(channel ==0)
			dst[index] = src[index].x;
		else if(channel ==1)
			dst[index] = src[index].y;
		else if(channel ==2)
			dst[index] = src[index].z;
	}
}

void CopyFloatFromFloat4(float4 *src, float * dst, const unsigned int width,const unsigned int height,const unsigned int channel)
{
	unsigned int  blocksW = (unsigned int) ceilf( (float) width / (float) BLOCK_SIZE );
    unsigned int  blocksH = (unsigned int) ceilf( (float) height /(float) BLOCK_SIZE );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );

	CopyFloatFromFloat4_Kernel<<< gridDim, blockDim >>>(src,dst, width, height,channel);
}

__global__ void CopyFloatToFloat4_Kernel(float *src, float4 * dst, const unsigned int width,const unsigned int height, int channel)
{
	int global_id_x= blockIdx.x*BLOCK_SIZE + threadIdx.x;
	int global_id_y= blockIdx.y*BLOCK_SIZE + threadIdx.y;

	if( global_id_x < width && global_id_y < height ){
		unsigned int index = global_id_x + global_id_y*width;
		if(channel ==0)
			dst[index].x = src[index];
		else if(channel ==1)
			dst[index].y = src[index];
		else if(channel ==2)
			dst[index].z = src[index];
	}
}

void CopyFloatToFloat4(float *src, float4 * dst, const unsigned int width,const unsigned int height,const unsigned int channel)
{
	unsigned int  blocksW = (unsigned int) ceilf( (float) width / (float) BLOCK_SIZE );
    unsigned int  blocksH = (unsigned int) ceilf( (float) height /(float) BLOCK_SIZE );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );

	CopyFloatToFloat4_Kernel<<< gridDim, blockDim >>>(src,dst, width, height,channel);
}

struct sax_functor {
	const float a;
	sax_functor(float p_a) : a(p_a) {}
	__host__ __device__ float operator()(const float& x) const {
		return a * x; } 
};

__global__ void ScaleVectors_Kernel(float2 * vector_buffer, float * norm_buffer, const float scale_value, int width, int height){
	
	int global_id_x= blockIdx.x*BLOCK_SIZE + threadIdx.x;
	int global_id_y= blockIdx.y*BLOCK_SIZE + threadIdx.y;

	if( global_id_x < width &&  global_id_y < height ){
		int write_pos = global_id_x + global_id_y*width;
		vector_buffer[write_pos].x*=scale_value;
		vector_buffer[write_pos].y*=scale_value;
		norm_buffer[write_pos]*=scale_value;
	}
}

__global__ void GeneralScaleVectors_Kernel(float2 * vector_buffer, float * norm_buffer, int width, int height){
	
	int global_id_x= blockIdx.x*BLOCK_SIZE + threadIdx.x;
	int global_id_y= blockIdx.y*BLOCK_SIZE + threadIdx.y;

	if( global_id_x < width &&  global_id_y < height ){
		int write_pos = global_id_x + global_id_y*width;
		float2 vector_value = vector_buffer[write_pos];
		float initial_norm = vector_value.x*vector_value.x + vector_value.y*vector_value.y;
		if(initial_norm > 0.f){
			initial_norm = sqrt(initial_norm);
			float new_norm = norm_buffer[write_pos];
			new_norm/=initial_norm;
			vector_buffer[write_pos].x*=new_norm;
			vector_buffer[write_pos].y*=new_norm;
		}
	}
}

void GLTools::GeneralVectorScaling(float2 * vector_buffer, float * norm_buffer, int width, int height){
	unsigned int  blocksW = (unsigned int) ceilf( (float) width / (float) BLOCK_SIZE );
    unsigned int  blocksH = (unsigned int) ceilf( (float) height /(float) BLOCK_SIZE );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );

	GeneralScaleVectors_Kernel<<< gridDim, blockDim >>>(vector_buffer,norm_buffer,width,height);
}

void GLTools::NormalizeVectorFieldByMaxima(float2 * vector_buffer, float * norm_buffer, int width, int height){

	thrust::device_ptr<float> dev_ptr_norms(norm_buffer);
	float max_value = thrust::reduce(dev_ptr_norms, dev_ptr_norms + height*width, (float) -FLT_MAX, thrust::maximum<float>());

	unsigned int  blocksW = (unsigned int) ceilf( (float) width / (float) BLOCK_SIZE );
    unsigned int  blocksH = (unsigned int) ceilf( (float) height /(float) BLOCK_SIZE );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );

	ScaleVectors_Kernel<<< gridDim, blockDim >>>(vector_buffer,norm_buffer,1.f/max_value,width,height);
}

///////////////////////////////////////// CUDA ADVECTION OBJECT ///////////////////////////////////////// 
#include "array-advection.h"

void CUDA_Advection_Object::SampleFieldFromGlobal_CUDA( VectorFieldMode vector_mode,float pos_w, float pos_h, bool normalized_coordinates, const float scale, NormalizationMode normalization_mode)
{
	scalar_texture.addressMode[0] = cudaAddressModeMirror;
	scalar_texture.addressMode[1] = cudaAddressModeMirror;
	scalar_texture.filterMode = cudaFilterModeLinear;
	scalar_texture.normalized = false;

	if(normalized_coordinates){
		pos_w *= (float)(width);
		pos_h *= (float)(height);
	}

	cutilSafeCall( cudaBindTextureToArray(scalar_texture,global_scalar_field_array) );

	unsigned int  blocksW = (unsigned int) ceilf( (float) width / (float) BLOCK_SIZE );
    unsigned int  blocksH = (unsigned int) ceilf( (float) height /(float) BLOCK_SIZE );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );

	bool normalize_field = normalization_mode == UNIFORM_NORMALIZATION ? true : false;

	if(vector_mode == MEAN_GRADIENT_FIELD){
		SampleGlobalMeanGradientField_Kernel<<< gridDim, blockDim >>>(float2_buffer, float_buffer_0, width, height,  1.f/(float)width, 1.f/(float)height, pos_w, pos_h, scale,normalize_field);
	}
	else if(vector_mode == ADVECTION_VECTOR_FIELD){
		SampleGlobalAdvectionField_Kernel<<< gridDim, blockDim >>>( float2_buffer, float_buffer_0, width, height,  1.f/(float)width, 1.f/(float)height, pos_w, pos_h, scale,normalize_field);
	}

	if(normalization_mode == MAXIMA_NORMALIZATION){
		GLTools::NormalizeVectorFieldByMaxima(float2_buffer,float_buffer_0,width,height);
	}

    cutilSafeCall( cudaUnbindTexture( scalar_texture ) );
}

void CUDA_Advection_Object::SampleFieldFromLocal_CUDA( VectorFieldMode vector_mode, NormalizationMode normalization_mode)
{
	scalar_texture.addressMode[0] = cudaAddressModeMirror;
	scalar_texture.addressMode[1] = cudaAddressModeMirror;
	scalar_texture.filterMode = cudaFilterModeLinear;
	scalar_texture.normalized = false;

	cutilSafeCall( cudaBindTextureToArray(scalar_texture,local_scalar_field_array) );

	unsigned int  blocksW = (unsigned int) ceilf( (float) width / (float) BLOCK_SIZE );
    unsigned int  blocksH = (unsigned int) ceilf( (float) height /(float) BLOCK_SIZE );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );

	bool normalize_field = normalization_mode == UNIFORM_NORMALIZATION ? true : false;

	if(vector_mode == MEAN_GRADIENT_FIELD){
		SampleLocalMeanGradientField_Kernel<<< gridDim, blockDim >>>(float2_buffer, float_buffer_0, width,height,normalize_field);
	}
	else if(vector_mode == MEAN_GRADIENT_FIELD_LINEAR){
		SampleLocalMeanGradientFieldLinear_Kernel<<< gridDim, blockDim >>>(float2_buffer, float_buffer_0, width,height,normalize_field);
	}
	else if(vector_mode == ADVECTION_VECTOR_FIELD){
		SampleLocalAdvectionField_Kernel<<< gridDim, blockDim >>>(float2_buffer, float_buffer_0, width, height,normalize_field);
	}

	if(normalization_mode == MAXIMA_NORMALIZATION){
		GLTools::NormalizeVectorFieldByMaxima(float2_buffer,float_buffer_0,width,height);
	}

    cutilSafeCall( cudaUnbindTexture( scalar_texture ) );
}


void CUDA_Advection_Object::SampleField_CUDA(VectorFieldMode p_vector_field_mode, SamplingDomain p_sampling_domain_mode, float pos_w, float pos_h, bool normalized_coordinates, const float scale, NormalizationMode normalization_mode)
{
	if(p_vector_field_mode == MEAN_GRADIENT_FIELD || p_vector_field_mode == ADVECTION_VECTOR_FIELD){
		if(p_sampling_domain_mode == GLOBAL_DOMAIN){
		SampleFieldFromGlobal_CUDA(p_vector_field_mode,pos_w, pos_h, normalized_coordinates, scale, normalization_mode);
		}
		else if(p_sampling_domain_mode == LOCAL_DOMAIN){
		SampleFieldFromLocal_CUDA( p_vector_field_mode,normalization_mode);
		}
		else if(p_sampling_domain_mode == HYBRID_DOMAIN){
		//SampleFieldFromLocal_CUDA( MEAN_GRADIENT_FIELD,normalization_mode);
		SampleFieldFromLocal_CUDA( MEAN_GRADIENT_FIELD_LINEAR,normalization_mode);
		cudaMemcpy(float_buffer_1,float_buffer_0,width*height*sizeof(float),cudaMemcpyDeviceToDevice);
		//SampleFieldFromGlobal_CUDA(p_vector_field_mode,pos_w, pos_h, normalized_coordinates,scale,NONE_NORMALIZATION);
		//SampleGlobalFieldFromSampledFilter_CUDA(p_vector_field_mode,pos_w, pos_h, normalized_coordinates,scale,NONE_NORMALIZATION);
		SampleGlobalFieldFromConstantFilter_CUDA(p_vector_field_mode,pos_w, pos_h, normalized_coordinates,scale,NONE_NORMALIZATION);
		GLTools::GeneralVectorScaling(float2_buffer,float_buffer_1,width,height);
		}

	}
}

void CUDA_Advection_Object::UpdateAdvectionArray_CUDA(SamplingDomain p_sampling_domain_mode,float pos_w, float pos_h, bool normalized_coordinates, const float scale, NormalizationMode normalization_mode)
{
	SampleField_CUDA(ADVECTION_VECTOR_FIELD,p_sampling_domain_mode,pos_w, pos_h, normalized_coordinates, scale, normalization_mode);
	cutilSafeCall( cudaMemcpyToArray(advection_field_array, 0, 0, float2_buffer, width*height*sizeof(float2), cudaMemcpyDeviceToDevice ) );
}

void CUDA_Advection_Object::BackwardSampling_CUDA(SamplingDomain color_sampling_domain, cudaArray * p_color_field_array, cudaArray * p_advection_field_array, const unsigned int num_iter, const unsigned int supersampling, float pos_w, float pos_h, bool normalized_coordinates, const float scale, const float step_amplification)
{
	unsigned int  blocksW = (unsigned int) ceilf( (float) width / (float) BLOCK_SIZE );
    unsigned int  blocksH = (unsigned int) ceilf( (float) height /(float) BLOCK_SIZE );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );

	color_texture.addressMode[0] = cudaAddressModeMirror;
	color_texture.addressMode[1] = cudaAddressModeMirror;
	color_texture.filterMode = cudaFilterModeLinear;
	color_texture.normalized = false;

	cutilSafeCall( cudaBindTextureToArray(color_texture,p_color_field_array) );

	advection_texture.addressMode[0] = cudaAddressModeMirror;
	advection_texture.addressMode[1] = cudaAddressModeMirror;
	advection_texture.filterMode = cudaFilterModeLinear;
	advection_texture.normalized = false; // NOTE : This may changed if out of range coordinates are used!!

	cutilSafeCall( cudaBindTextureToArray(advection_texture,p_advection_field_array) );

	if(normalized_coordinates){
		pos_w *= (float)(width);
		pos_h *= (float)(height);
	}

	//if(level_number==2){
	//	pos_w = 0.0000f;
	//	pos_h = 0.0000f;
	//	scale = 1.f;
	//}

	if(color_sampling_domain == GLOBAL_DOMAIN)
	BackwardSamplingFromGlobal_Kernel<<< gridDim, blockDim >>>(float4_buffer, width, height, pos_w, pos_h, scale, num_iter, supersampling, 0.f,step_amplification);
	else if(color_sampling_domain == LOCAL_DOMAIN)
	BackwardSamplingFromLocal_Kernel<<< gridDim, blockDim >>>(float4_buffer, width, height, num_iter, supersampling, 0.f,step_amplification);

	cutilSafeCall( cudaUnbindTexture(color_texture) );
    cutilSafeCall( cudaUnbindTexture(advection_texture) );
}


__global__ void SampleGlobalAdvectionFromSampledFilter_Kernel( float2* dst, float * norms, const unsigned int imgWidth,const unsigned int imgHeight, const float finv_imgWidht,  const float finv_imgHeight, const float corner_w,const float corner_h,const float scale,  bool p_normalize_advection_field, float * filter_w, float * filter_h, const int filter_support, const int filter_radius)
{
	int thread_idx = blockIdx.x*BLOCK_SIZE + threadIdx.x;
	int thread_idy = blockIdx.y*BLOCK_SIZE + threadIdx.y;

	if( thread_idx < imgWidth && thread_idy < imgHeight ){

	float dw = corner_w + scale*((float)thread_idx + 0.5f);//global_pos_w
	float dh = corner_h + scale*((float)thread_idy + 0.5f);//global_pos_h
	
	float cw = floor(dw);
	float ch = floor(dh);

	dw-=cw;
	dh-=ch;

	float filter_pos_w;
	float filter_pos_h;

	float grad_x = 0.f; 
	float grad_y = 0.f;
	float hessian_xx = 0.f;
	float hessian_yy = 0.f;
	float hessian_xy = 0.f;

	float filter_value_h;
	float filter_value;
	float scalar_value;

	float shift_h;

	for (int ih = 0; ih < filter_support; ih++){
		filter_value_h = filter_h[imgHeight*ih + thread_idy];//texFetch may be more efficient!!
			shift_h = (float)(ih - filter_radius + 1);
			for(int iw = 0 ; iw < filter_support; iw++){
				//texFetch may be more efficient!!
				filter_pos_w = dw - (float)(iw - filter_radius + 1); // CONVOLUTION POSITION // shift_w
				filter_pos_h = dh - shift_h;
				scalar_value = tex2D(scalar_texture,cw + (float)(iw - filter_radius + 1),ch + shift_h);
				filter_value = filter_value_h*filter_w[imgWidth*iw + thread_idx];//texFetch may be more efficient!!

				grad_x +=(-filter_pos_w / SIGMA_SQUARED)*filter_value*scalar_value;
				grad_y +=(-filter_pos_h / SIGMA_SQUARED)*filter_value*scalar_value;
				hessian_xx += ((filter_pos_w*filter_pos_w - SIGMA_SQUARED)/(SIGMA_SQUARED*SIGMA_SQUARED))*filter_value*scalar_value;
				hessian_yy += ((filter_pos_h*filter_pos_h - SIGMA_SQUARED)/(SIGMA_SQUARED*SIGMA_SQUARED))*filter_value*scalar_value;
				hessian_xy += ((filter_pos_w*filter_pos_h)/(SIGMA_SQUARED*SIGMA_SQUARED))*filter_value*scalar_value;
			}
	}

	scalar_value =hessian_xx*grad_x*grad_x + 2.f*hessian_xy*grad_x*grad_y +hessian_yy*grad_y*grad_y;
	grad_x *= scalar_value;
	grad_y *= scalar_value;

	float grad_norm = sqrt(grad_x*grad_x + grad_y*grad_y);

	if(p_normalize_advection_field && grad_norm > 0.f){
			grad_x/=grad_norm;
			grad_y/=grad_norm;
			norms[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth] = 1.f;
	}
	else{
		norms[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth] = grad_norm;
	}

	dst[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth].x  = grad_x;
	dst[(blockIdx.x*BLOCK_SIZE + threadIdx.x) + (blockIdx.y*BLOCK_SIZE + threadIdx.y)*imgWidth].y  = grad_y;
	}
}

void CUDA_Advection_Object::SampleGlobalFieldFromSampledFilter_CUDA(VectorFieldMode vector_mode,float pos_w, float pos_h, bool normalized_coordinates, const float scale, NormalizationMode normalization_mode)
{
	scalar_texture.addressMode[0] = cudaAddressModeMirror;
	scalar_texture.addressMode[1] = cudaAddressModeMirror;
	scalar_texture.filterMode = cudaFilterModeLinear;
	scalar_texture.normalized = false;

	if(normalized_coordinates){
		pos_w *= (float)(width);
		pos_h *= (float)(height);
	}

	cutilSafeCall( cudaBindTextureToArray(scalar_texture,global_scalar_field_array) );

	unsigned int  blocksW = (unsigned int) ceilf( (float) width / (float) BLOCK_SIZE );
    unsigned int  blocksH = (unsigned int) ceilf( (float) height /(float) BLOCK_SIZE );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );

	bool normalize_field = normalization_mode == UNIFORM_NORMALIZATION ? true : false;

	if(vector_mode == MEAN_GRADIENT_FIELD){
		printf("unimplemented!! \n");		
	}
	else if(vector_mode == ADVECTION_VECTOR_FIELD){
		//SampleGlobalAdvectionField_Kernel<<< gridDim, blockDim >>>( float2_buffer, float_buffer_0, width, height,  1.f/(float)width, 1.f/(float)height, pos_w, pos_h, scale,normalize_field);
	    SampleGlobalAdvectionFromSampledFilter_Kernel<<< gridDim, blockDim >>>(float2_buffer, float_buffer_0, width, height,  1.f/(float)width, 1.f/(float)height, pos_w, pos_h, scale,normalize_field,sampled_filter_values_w,sampled_filter_values_h,filter_support,filter_radius);
	}
	if(normalization_mode == MAXIMA_NORMALIZATION){
		GLTools::NormalizeVectorFieldByMaxima(float2_buffer,float_buffer_0,width,height);
	}
    cutilSafeCall( cudaUnbindTexture( scalar_texture ) );
}

void CUDA_Advection_Object::SampleGlobalFieldFromConstantFilter_CUDA(VectorFieldMode vector_mode,float pos_w, float pos_h, bool normalized_coordinates, const float scale, NormalizationMode normalization_mode)
{
	scalar_texture.addressMode[0] = cudaAddressModeMirror;
	scalar_texture.addressMode[1] = cudaAddressModeMirror;
	scalar_texture.filterMode = cudaFilterModeLinear;
	scalar_texture.normalized = false;

	if(normalized_coordinates){
		pos_w *= (float)(width);
		pos_h *= (float)(height);
	}

	cutilSafeCall( cudaBindTextureToArray(scalar_texture,global_scalar_field_array) );

	unsigned int  blocksW = (unsigned int) ceilf( (float) width / (float) BLOCK_SIZE );
    unsigned int  blocksH = (unsigned int) ceilf( (float) height /(float) BLOCK_SIZE );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );

	bool normalize_field = normalization_mode == UNIFORM_NORMALIZATION ? true : false;

	if(vector_mode == MEAN_GRADIENT_FIELD){
		printf("unimplemented!! \n");		
	}
	else if(vector_mode == ADVECTION_VECTOR_FIELD){
		//SampleGlobalAdvectionField_Kernel<<< gridDim, blockDim >>>( float2_buffer, float_buffer_0, width, height,  1.f/(float)width, 1.f/(float)height, pos_w, pos_h, scale,normalize_field);
		SampleGlobalAdvectionField_Constant_Kernel<<< gridDim, blockDim >>>(float2_buffer, float_buffer_0, width, height,  1.f/(float)width, 1.f/(float)height, pos_w, pos_h, scale,normalize_field);
	}
	if(normalization_mode == MAXIMA_NORMALIZATION){
		GLTools::NormalizeVectorFieldByMaxima(float2_buffer,float_buffer_0,width,height);
	}
    cutilSafeCall( cudaUnbindTexture( scalar_texture ) );
}

__global__ void FilterSampling_Kernel( float * sampled_values_w,float * sampled_values_h,const int width, const int height, const int  filter_support, const int  filter_radius, const float corner_w, const float corner_h, const float scale)
{
	int thread_id = blockIdx.x*LARGE_BLOCK_SIZE + threadIdx.x;
	float shift;

	if( thread_id < width ){
		shift = corner_w+ scale*((float)thread_id + 0.5f); //global_pos
		shift -= floor(shift);
		shift += (filter_radius - 1.f);
		for(int i=0; i<filter_support; i++){
			sampled_values_w[width*i + thread_id] = Gaussian(shift);
			shift -=1.f;
		}
	}

	if( thread_id < height ){
		shift = corner_h+ scale*((float)thread_id + 0.5f); //global_pos
		shift -= floor(shift);
		shift += (filter_radius - 1.f);
		for(int i=0; i<filter_support; i++){
			sampled_values_h[height*i + thread_id] = Gaussian(shift);
			shift -=1.f;
		}
	}
}

void  CUDA_Advection_Object::UpdateFilterSamples(float pos_w, float pos_h, bool normalized_coordinates, const float scale)
{
	int largest_dimension = width > height ? width : height;
	unsigned int  blocks = (unsigned int) ceilf( (float) largest_dimension / (float) LARGE_BLOCK_SIZE );
    dim3 gridDim( blocks, 1, 1 );
    dim3 blockDim( LARGE_BLOCK_SIZE, 1, 1 );

	if(normalized_coordinates){
		pos_w *= (float)(width);
		pos_h *= (float)(height);
	}

	FilterSampling_Kernel<<< gridDim, blockDim >>>(sampled_filter_values_w,sampled_filter_values_h,width, height, filter_support, filter_radius, pos_w, pos_h, scale);
}


///////////////////////////////////////// ARRAY ADVECTION  ///////////////////////////////////////// 

///////////////////////////////////////// ARRAY HIERARCHY ///////////////////////////////////////// 

#include "array-hierarchy.h"

void Array_Hierarchy::SetTextureToLevel(const int level_number,cudaGraphicsResource_t& p_texture_resource)
{
	cutilSafeCall( cudaMemcpyFromArray(float4_buffer,hierarchical_color_field_array[level_number],0, 0, width*height*sizeof(float4), cudaMemcpyDeviceToDevice ) );
	GLTools::TransformFloat4ArrayToGLColorBuffer_CUDA(p_texture_resource,float4_buffer,uchar4_buffer,width,height);
}

void Array_Hierarchy::SamplePreviousLevel_CUDA(const int level_number){

	// Sample Color and Scalar Field
	color_texture.addressMode[0] = cudaAddressModeMirror;
	color_texture.addressMode[1] = cudaAddressModeMirror;
	color_texture.filterMode = cudaFilterModeLinear;
	color_texture.normalized = false;

	cutilSafeCall( cudaBindTextureToArray(color_texture,hierarchical_color_field_array[level_number-1]));

	unsigned int  blocksW = (unsigned int) ceilf( (float) width / (float) BLOCK_SIZE );
    unsigned int  blocksH = (unsigned int) ceilf( (float) height /(float) BLOCK_SIZE );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );

	Point<2> relative_position = (level_absolute_coordinate[level_number] - level_absolute_coordinate[level_number-1])/level_absolute_scale[level_number-1];
	float pos_w = (float)relative_position[0]*(float)(width-1);
	float pos_h = (float)relative_position[1]*(float)(height-1);
	float relative_scale = level_absolute_scale[level_number]/level_absolute_scale[level_number-1];

	SamplePreviousLevel_Kernel<<< gridDim, blockDim >>>(float4_buffer,float_buffer_0,pos_w,pos_h,relative_scale,width,height);

    cutilSafeCall( cudaUnbindTexture( color_texture ) );
}

void Array_Hierarchy::UpdateLevel(const int num_iter, const int level_number)
{
	SamplePreviousLevel_CUDA(level_number);
	cutilSafeCall( cudaMemcpyToArray(local_color_field_array, 0, 0, float4_buffer, width*height*sizeof(float4), cudaMemcpyDeviceToDevice ) );
	cutilSafeCall( cudaMemcpyToArray(local_scalar_field_array, 0, 0, float_buffer_0, width*height*sizeof(float), cudaMemcpyDeviceToDevice ) );
	
	UpdateAdvectionArray_CUDA(HYBRID_DOMAIN,(float)level_absolute_coordinate[level_number][0],(float)level_absolute_coordinate[level_number][1],true,(float)level_absolute_scale[level_number],MAXIMA_NORMALIZATION);
	float step_advection = 1.f;
	
	BackwardSampling_CUDA(LOCAL_DOMAIN,local_color_field_array,advection_field_array,10,3,0.f,0.f,true,0.f,step_advection);
	cutilSafeCall( cudaMemcpyToArray(hierarchical_color_field_array[level_number], 0, 0, float4_buffer, width*height*sizeof(float4), cudaMemcpyDeviceToDevice ) );
}

////////////////////////////////////////////// GRADIENT ACCUMULATION
__global__ void AccumulateGradient_Kernel(float3 * dst_grad_w,float3 * dst_grad_h,const unsigned int width,const unsigned int height, const unsigned int iter_num,const unsigned int subdiv,const unsigned int it_w,const unsigned int it_h, const float step_amplification)
{

	int index_w = (blockIdx.x*subdiv + it_w)*BLOCK_SIZE_X + threadIdx.x;
 	int index_h = (blockIdx.y*subdiv + it_h)*BLOCK_SIZE_Y + threadIdx.y;

	if( (index_w) < width && (index_h) < height ){

        float dual_pos_w = (float)(index_w) + 0.5f;
		float dual_pos_h = (float)(index_h) + 0.5f;
		float4 center = tex2D(color_texture,dual_pos_w,dual_pos_h); 
		float2 advection_direction;

		// Advect grad_w
		if( index_w < width-1){
			
			float pos_w = dual_pos_w + 0.5f;
			float pos_h = dual_pos_h;

			for(int i =0; i<iter_num; i++){
				advection_direction = tex2D(advection_texture,pos_w,pos_h);
				pos_w += advection_direction.x*step_amplification;
				pos_h += advection_direction.y*step_amplification;
			}

			if((pos_w >= 1.f && pos_w <= (float)(width) - 1.f) && (pos_h >= 0.5f && pos_h <= (float)(height) -0.5f)){
				int cw = (int)floor(pos_w);
				int ch = (int)floor(pos_h-0.5f);
				float dw = (pos_w)-(float)(cw);
				float dh = (pos_h-0.5f)-(float)(ch);

				float4 neighbour = tex2D(color_texture,dual_pos_w + 1.f,dual_pos_h); 
				float3 grad = SubstractionFloat4ToFloat3(neighbour,center);
				for(int i=0; i<BLOCK_SIZE_X; i++){
					for(int j=0; j<BLOCK_SIZE_Y; j++){
					if(threadIdx.x == i && threadIdx.y ==j){

						dst_grad_w[cw+ch*(width+1)].x += grad.x*(1.f-dw)*(1.f-dh);
						dst_grad_w[cw+ch*(width+1)].y += grad.y*(1.f-dw)*(1.f-dh);
						dst_grad_w[cw+ch*(width+1)].z += grad.z*(1.f-dw)*(1.f-dh);

						dst_grad_w[cw+1+ch*(width+1)].x += grad.x*dw*(1.f-dh);
						dst_grad_w[cw+1+ch*(width+1)].y += grad.y*dw*(1.f-dh);
						dst_grad_w[cw+1+ch*(width+1)].z += grad.z*dw*(1.f-dh);

						dst_grad_w[cw + (ch+1)*(width+1)].x += grad.x*(1.f-dw)*dh;
						dst_grad_w[cw + (ch+1)*(width+1)].y += grad.y*(1.f-dw)*dh;
						dst_grad_w[cw + (ch+1)*(width+1)].z += grad.z*(1.f-dw)*dh;

						dst_grad_w[cw + 1 + (ch+1)*(width+1)].x += grad.x*dw*dh;
						dst_grad_w[cw + 1 + (ch+1)*(width+1)].y += grad.y*dw*dh;
						dst_grad_w[cw + 1 + (ch+1)*(width+1)].z += grad.z*dw*dh;
						}
					}
				}
			}
		}

		if(index_h < height-1){
		
			float pos_w = dual_pos_w;
			float pos_h = dual_pos_h + 0.5f;

			for(int i =0; i<iter_num; i++){
				advection_direction = tex2D(advection_texture,pos_w,pos_h);
				pos_w += advection_direction.x*step_amplification;
				pos_h += advection_direction.y*step_amplification;
			}

			if((pos_w >= 0.5 && pos_w <= (float)(width) -0.5f) && (pos_h >= 1.f && pos_h <= (float)(height) -1.f)){
				int cw = (int)floor(pos_w-0.5f);
				int ch = (int)floor(pos_h);
				float dw = (pos_w)-(float)(cw-0.5f);
				float dh = (pos_h)-(float)(ch);

				float4 neighbour = tex2D(color_texture,dual_pos_w,dual_pos_h + 1.f); 
				float3 grad = SubstractionFloat4ToFloat3(neighbour,center);
				for(int i=0; i<BLOCK_SIZE_X; i++){
					for(int j=0; j<BLOCK_SIZE_Y; j++){
					if(threadIdx.x == i && threadIdx.y ==j){
						dst_grad_h[cw+ch*(width)].x += grad.x*(1.f-dw)*(1.f-dh);
						dst_grad_h[cw+ch*(width)].y += grad.y*(1.f-dw)*(1.f-dh);
						dst_grad_h[cw+ch*(width)].z += grad.z*(1.f-dw)*(1.f-dh);

						dst_grad_h[cw+1+ch*(width)].x += grad.x*dw*(1.f-dh);
						dst_grad_h[cw+1+ch*(width)].y += grad.y*dw*(1.f-dh);
						dst_grad_h[cw+1+ch*(width)].z += grad.z*dw*(1.f-dh);


						dst_grad_h[cw + (ch+1)*(width)].x += grad.x*(1.f-dw)*dh;
						dst_grad_h[cw + (ch+1)*(width)].y += grad.y*(1.f-dw)*dh;
						dst_grad_h[cw + (ch+1)*(width)].z += grad.z*(1.f-dw)*dh;

						dst_grad_h[cw + 1 + (ch+1)*(width)].x += grad.x*dw*dh;
						dst_grad_h[cw + 1 + (ch+1)*(width)].y += grad.y*dw*dh;
						dst_grad_h[cw + 1 + (ch+1)*(width)].z += grad.z*dw*dh;
						}
					}
				}
			}
		}
	}
}

__global__ void Divergence_Kernel(float3 * src_w,float3 * src_h,float* dst, const unsigned int imgWidth,const unsigned int imgHeight, int channel)// COMPUTES MINUS DIVERGENCE
{
	int global_id_x= blockIdx.x*BLOCK_SIZE + threadIdx.x;
	int global_id_y= blockIdx.y*BLOCK_SIZE + threadIdx.y;
	if( global_id_x < imgWidth && global_id_y < imgHeight){
			int write_pos = global_id_x + global_id_y*imgWidth;
			if(channel ==0)
				dst[write_pos] =  (src_w[write_pos+global_id_y].x - src_w[write_pos+global_id_y+1].x) + (src_h[write_pos].x - src_h[write_pos + imgWidth].x);
			else if(channel == 1)
				dst[write_pos] =  (src_w[write_pos+global_id_y].y - src_w[write_pos+global_id_y+1].y) + (src_h[write_pos].y - src_h[write_pos + imgWidth].y);
			else if(channel == 2)
				dst[write_pos] =  (src_w[write_pos+global_id_y].z - src_w[write_pos+global_id_y+1].z) + (src_h[write_pos].z - src_h[write_pos + imgWidth].z);
	}
}

void CUDA_Advection_Object::Divergence_CUDA(float3 * src_w,float3 * src_h,float* dst, const unsigned int imgWidth,const unsigned int imgHeight, int channel)// COMPUTES MINUS DIVERGENCE
{
	unsigned int  blocksW = (unsigned int) ceilf( (float) width / (float) BLOCK_SIZE );
    unsigned int  blocksH = (unsigned int) ceilf( (float) height /(float) BLOCK_SIZE );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );

	Divergence_Kernel<<< gridDim, blockDim >>>(gradient_accumulation_array_w,gradient_accumulation_array_h,float_buffer_0,width,height,channel);
}

void CUDA_Advection_Object::AccumulateGradient_CUDA(const unsigned int num_iter , const float step_amplification)
{

	cutilSafeCall(cudaMemset(gradient_accumulation_array_w, 0, ( width+1 )* height * sizeof(float3)));
	cutilSafeCall(cudaMemset(gradient_accumulation_array_h, 0, width * (height +1) * sizeof(float3)));

	color_texture.addressMode[0] = cudaAddressModeMirror;
	color_texture.addressMode[1] = cudaAddressModeMirror;
	color_texture.filterMode = cudaFilterModeLinear;
	color_texture.normalized = false;

	cutilSafeCall( cudaBindTextureToArray(color_texture,local_color_field_array) );

	advection_texture.addressMode[0] = cudaAddressModeMirror;
	advection_texture.addressMode[1] = cudaAddressModeMirror;
	advection_texture.filterMode = cudaFilterModeLinear;
	advection_texture.normalized = false; // NOTE : This may changed if out of range coordinates are used!!

	cutilSafeCall( cudaBindTextureToArray(advection_texture,advection_field_array) );

	const unsigned int subdiv =2;

	const unsigned int num_threads_w = ceilf( (float) width / (float) (subdiv) );
	const unsigned int num_threads_h = ceilf( (float) height / (float) (subdiv) );

	const unsigned int  blocksW = (unsigned int) ceilf( (float) num_threads_w / (float) BLOCK_SIZE_X );
    const unsigned int  blocksH = (unsigned int) ceilf( (float) num_threads_h /(float) BLOCK_SIZE_Y );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE_X, BLOCK_SIZE_Y, 1 );

	for(int it_w =0; it_w<subdiv; it_w++){
			for(int it_h =0; it_h<subdiv; it_h++){
				AccumulateGradient_Kernel<<< gridDim, blockDim >>>(gradient_accumulation_array_w,gradient_accumulation_array_h, width, height,num_iter,subdiv,it_w,it_h,step_amplification);
			}
	}
	cutilSafeCall( cudaUnbindTexture(color_texture) );
    cutilSafeCall( cudaUnbindTexture(advection_texture) );
		
}


__global__ void Gradient_Kernel(float3 * gradient_accumulation_array_w, float3 * gradient_accumulation_array_h, const int width, const int height)
{
	int global_id_x= blockIdx.x*BLOCK_SIZE + threadIdx.x;
	int global_id_y= blockIdx.y*BLOCK_SIZE + threadIdx.y;
	if( global_id_x < width  && global_id_y < height){

		float pos_w = (float)(global_id_x) + 0.5f;
		float pos_h = (float)(global_id_y) + 0.5f;

		int write_pos  = global_id_x + global_id_y*width;
		float4 center = tex2D(color_texture,pos_w,pos_h); 

		if( global_id_x < width-1){
		
		float4 neighbour_right = tex2D(color_texture,pos_w + 1.f,pos_h); 
		write_pos = (global_id_x + 1) + global_id_y*(width+1);

		gradient_accumulation_array_w[write_pos] = SubstractionFloat4ToFloat3(neighbour_right,center);
		}

		if(global_id_y < height-1){
		
		float4 neighbour_up = tex2D(color_texture,pos_w,pos_h + 1.f); 
		write_pos = global_id_x+ (global_id_y + 1)*width;
		gradient_accumulation_array_h[write_pos] = SubstractionFloat4ToFloat3(neighbour_up,center);
		}
	}
}

void CUDA_Advection_Object::SolveGradientField_CUDA()
{
	for(int i=0; i<3; i++){
		 CopyFloatFromFloat4(float4_buffer,float_buffer_0,width,height,i);
		 thrust::device_ptr<float> dev_ptr(float_buffer_0);
		 float dc = thrust::reduce(dev_ptr, dev_ptr + width * height, (float) 0.f ,thrust::plus<float>()); 
		 Divergence_CUDA(gradient_accumulation_array_w,gradient_accumulation_array_h,float_buffer_0,width,height,i);
		 GPUSolvers::FFTLaplaceSolver(float_buffer_0,fft_padding_buffer,fft_complex_buffer,fftPlanFwd, fftPlanInv,width,height,dc);
		 CopyFloatToFloat4(float_buffer_0,float4_buffer,width,height,i);
	}
}

