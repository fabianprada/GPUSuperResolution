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


#ifndef CUDA_ADVECTION_INCLUDED
#define CUDA_ADVECTION_INCLUDED
#include "zoom-visualization.h"
#include <cufft.h>

extern "C" void SetGaussianInteger();

class CUDA_Advection_Object : public ZoomInObject{
public:
	CUDA_Advection_Object(int p_height, int p_width, unsigned char * p_color, Canvas_Rectangle * p_canvas){
		color = p_color;
		height = p_height;
		width = p_width;
		canvas = p_canvas;
		vector_field_mode =  NONE_VECTOR_FIELD;
		filter_radius = 6;
		filter_support = 2*filter_radius;
	}

	// Image Handle
	unsigned char * color;
	int height;
	int width;
	Canvas_Rectangle *  canvas;

	// GL Buffers
	GLuint texture_buffer;
	GLuint vector_field_buffer;

	// Interop Resources
	cudaGraphicsResource_t texture_resource;
	cudaGraphicsResource_t vector_field_resource;

	// CUDA Arrays
	cudaArray * color_field_array;
	cudaArray * global_scalar_field_array;
	cudaArray * local_scalar_field_array;
	cudaArray * local_color_field_array;
	cudaArray * advection_field_array;

	// Auxiliar CUDA Buffer
	float4 * float4_buffer;
	float2 * float2_buffer;
	float  *  float_buffer_0;
	float  *  float_buffer_1;

	unsigned char * uchar4_buffer;

	float3 * gradient_accumulation_array_w;
	float3 * gradient_accumulation_array_h;

	float  * fft_padding_buffer;
	float2  * fft_complex_buffer;

	cufftHandle fftPlanFwd;
	cufftHandle fftPlanInv;

	int filter_radius;
	int filter_support;
	float * sampled_filter_values_w;
	float * sampled_filter_values_h;

	//int advection_steps;
	Point<2> roi_absolute_coordinate;
	double roi_absolute_scale;

	void Initialize_CUDA_Advection_Object();
	VectorFieldMode vector_field_mode;

	void SampleFieldFromGlobal_CUDA( VectorFieldMode vector_mode, float pos_w, float pos_h, bool normalized_coordinates, const float scale, NormalizationMode normalization_mode);
	void SampleFieldFromLocal_CUDA( VectorFieldMode vector_mode, NormalizationMode normalization_mode);

	void SampleField_CUDA(VectorFieldMode p_vector_field_mode, SamplingDomain p_sampling_domain_mode,float pos_w, float pos_h, bool normalized_coordinates, const float scale, NormalizationMode normalization_mode);
	void UpdateAdvectionArray_CUDA(SamplingDomain p_sampling_domain_mode,float pos_w, float pos_h, bool normalized_coordinates, const float scale, NormalizationMode normalization_mode);

	void BackwardSampling_CUDA(SamplingDomain color_sampling_domain, cudaArray * p_color_field_array, cudaArray * p_advection_field_array, const unsigned int num_iter, const unsigned int supersampling, float pos_w, float pos_h, bool normalized_coordinates, const float scale, const float step_amplification);
	void ChangeVectorField_CallBack(Prompt  &  prompt);


	void SampleGlobalFieldFromSampledFilter_CUDA(VectorFieldMode vector_mode,float pos_w, float pos_h, bool normalized_coordinates, const float scale, NormalizationMode normalization_mode);
	void UpdateFilterSamples(float pos_w, float pos_h, bool normalized_coordinates, const float scale);

	void Divergence_CUDA(float3 * src_w,float3 * src_h,float* dst, const unsigned int imgWidth,const unsigned int imgHeight, int channel);// COMPUTES MINUS DIVERGENCE
	void AccumulateGradient_CUDA(const unsigned int num_iter , const float step_amplification);
	void SolveGradientField_CUDA();

	void SampleGlobalFieldFromConstantFilter_CUDA(VectorFieldMode vector_mode,float pos_w, float pos_h, bool normalized_coordinates, const float scale, NormalizationMode normalization_mode);



};

#endif //CUDA_ADVECTION_INCLUDED