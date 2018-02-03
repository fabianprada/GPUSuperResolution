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


#include "cuda-advection-object.h"
#include "cutil.h"
#include "cutil_inline_runtime.h"

void CUDA_Advection_Object::ChangeVectorField_CallBack(Prompt  &  prompt){
	if (vector_field_mode == NONE_VECTOR_FIELD){
		vector_field_mode = MEAN_GRADIENT_FIELD;
		sprintf(prompt.string_value, "VECTOR FIELD MODE : MEAN_GRADIENT_FIELD");
	}
	else if (vector_field_mode == MEAN_GRADIENT_FIELD){
		vector_field_mode = ADVECTION_VECTOR_FIELD;
		sprintf(prompt.string_value, "VECTOR FIELD MODE : ADVECTION_VECTOR_FIELD");
	}
	else if (vector_field_mode == ADVECTION_VECTOR_FIELD){
		vector_field_mode = NONE_VECTOR_FIELD;
		sprintf(prompt.string_value, "VECTOR FIELD MODE : NONE_VECTOR_FIELD");
	}
	else{
		sprintf(prompt.string_value, "Undefined visualization mode");
	}
	//vector_field_buffer_updated = false;
}

void CUDA_Advection_Object::Initialize_CUDA_Advection_Object()
{
	// Initialize GL Buffers
	GLTools::CreateVBO(vector_field_buffer,width * height * 9 * sizeof(float));
	GLTools::CreateTexture(texture_buffer,width,height);

	// Initialize CUDA Resources
	GLTools::CreateCUDABufferResource( vector_field_resource, vector_field_buffer, cudaGraphicsMapFlagsWriteDiscard );
	GLTools::CreateCUDAImageResource( texture_resource, texture_buffer, cudaGraphicsMapFlagsWriteDiscard );


	// Initialize CUDA Arrays
	cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindFloat); 
	cudaMallocArray(&global_scalar_field_array,&channelDesc1, width, height);

	cudaChannelFormatDesc channelDesc3 = cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindFloat); 
	cudaMallocArray(&local_scalar_field_array,&channelDesc3, width, height);
	
	cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc(32,32,32,32,cudaChannelFormatKindFloat); 
	cudaMallocArray(&color_field_array,&channelDesc2,width,height);
	cudaMallocArray(&local_color_field_array,&channelDesc2,width,height);

	cudaChannelFormatDesc advection_field_desc = cudaCreateChannelDesc(32,32,0,0,cudaChannelFormatKindFloat); 
	cudaMallocArray(&advection_field_array,&advection_field_desc, width, height);


	// Initialize CUDA Buffers
	cudaMalloc((void **)&float_buffer_0,width* height * sizeof(float));
	cudaMalloc((void **)&float_buffer_1,width* height * sizeof(float));
	cudaMalloc((void **)&float2_buffer,width* height * sizeof(float2));
	cudaMalloc((void **)&float4_buffer,width* height * sizeof(float4));
	cudaMalloc((void**)&uchar4_buffer,height*width*4*sizeof(char));

	cudaMalloc((void **)&sampled_filter_values_w,width*(2*filter_radius)* sizeof(float));
	cudaMalloc((void **)&sampled_filter_values_h,height*(2*filter_radius)* sizeof(float));

	// Initialize Values
	float * luminance = CastingTools::TransformUnsignedToNormalizedFloat(color,width*height,4,1);
	cudaMemcpyToArray(global_scalar_field_array, 0, 0, luminance,width*height*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpyToArray(local_scalar_field_array, 0, 0, luminance,width*height*sizeof(float), cudaMemcpyHostToDevice);
	delete luminance;

	float * color_values = CastingTools::TransformUnsignedToNormalizedFloat(color,width*height,4,4);
	cudaMemcpyToArray(color_field_array, 0, 0, color_values,width*height*4*sizeof(float), cudaMemcpyHostToDevice);
	delete color_values;

	unsigned int gradient_accumulation_w_size = (width+1) * height * sizeof(float3);
	cudaMalloc((void **)&gradient_accumulation_array_w,gradient_accumulation_w_size);

	unsigned int gradient_accumulation_h_size = (width) * (height+1) * sizeof(float3);
	cudaMalloc((void **)&gradient_accumulation_array_h,gradient_accumulation_h_size);

	unsigned int padded_array_width = 2*width;
	unsigned int padded_array_height = 2*height;

	cudaMalloc((void**)&fft_padding_buffer,padded_array_width*padded_array_height*sizeof(float));
	cudaMalloc((void**)&fft_complex_buffer,(padded_array_height )*(padded_array_width/ 2 + 1) *sizeof(float2));

	cufftSafeCall( cufftPlan2d(&fftPlanFwd, padded_array_height, padded_array_width, CUFFT_R2C) );
    cufftSafeCall( cufftPlan2d(&fftPlanInv, padded_array_height, padded_array_width, CUFFT_C2R) );

	SetGaussianInteger();
}