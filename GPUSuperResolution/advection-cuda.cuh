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


#ifndef ADVECTION_CUDA_INCLUDED
#define ADVECTION_CUDA_INCLUDED

#include <cuda_runtime_api.h>
#include "cutil.h"
#include "cutil_inline_runtime.h"
#include <thrust\device_ptr.h>
#include <thrust\transform_reduce.h>
#include <thrust\functional.h>
#include "gpu-solver.h"

#define LARGE_BLOCK_SIZE 128
#define BLOCK_SIZE 16
#define BLOCK_SIZE_X 8
#define BLOCK_SIZE_Y 4

texture<float, cudaTextureType2D,cudaReadModeElementType> scalar_texture;// try float4
texture<float4, cudaTextureType2D,cudaReadModeElementType> color_texture;// try float4
texture<float2, cudaTextureType2D,cudaReadModeElementType> gradient_texture;
texture<float2, cudaTextureType2D,cudaReadModeElementType> vector_field_texture;
texture<float2, cudaTextureType2D,cudaReadModeElementType> advection_texture;

__global__ void SampleAdvectionField_Kernel( float2* dst,const unsigned int imgWidth,const unsigned int imgHeight, const float finv_imgWidht,  const float finv_imgHeight, const float corner_w,const float corner_h,const float scale, const unsigned int supersampling, bool p_normalize_advection_field );
__device__ void NormalFloat2UChar(const float x, unsigned char & u);
__device__ unsigned char NormalFloat2UChar(const float x);


__device__  void AdditionFloat4 (float4 & in,const float4 a,const float4 b);
__device__  void SubstractionFloat4 (float4 & in,const float4 a,const float4 b);

#endif //ADVECTION_CUDA_INCLUDED