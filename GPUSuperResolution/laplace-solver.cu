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


#include "gpu-solver.h"
#include <cuda_runtime_api.h>
#include <cuda_runtime.h>
#include "cutil.h"
#include "cutil_inline_runtime.h"
#include <stdio.h>
#include "file-io.h"

#define BLOCK_SIZE 16
#define PI2 6.2831853

//typedef float2 fcomplex; 
//float * rhs_h = 0;
//float * rhs_d = 0;
//float * symmetric_extended_rhs_d = 0;
//float * symmetric_extended_rhs_h = 0;
//fcomplex * symmetric_extended_rhs_fft_h = 0;
//fcomplex * symmetric_extended_rhs_fft_d = 0;
//unsigned int array_size = 0;

__global__ void ExtendedSymmetricRHS(float *  rhs_d,float *  extended_rhs_d,const unsigned int imgWidth,const unsigned int imgHeight)
{
	unsigned int tx =blockIdx.x*BLOCK_SIZE + threadIdx.x;
	unsigned int ty =blockIdx.y*BLOCK_SIZE + threadIdx.y;
	if(  tx < imgWidth &&  ty < imgHeight ){
		float value = rhs_d[tx + imgWidth*ty];
		extended_rhs_d[tx + 2*imgWidth*ty] = value;
		extended_rhs_d[ 2*imgWidth - 1 - tx + 2*imgWidth*ty] = value;
		extended_rhs_d[tx + 2*imgWidth*(2*imgHeight - 1 - ty)] = value;
		extended_rhs_d[2*imgWidth - 1 - tx + 2*imgWidth*(2*imgHeight - 1 - ty)] = value;
	}
}

__global__ void ExtractRHS(float *  rhs_d,float *  extended_rhs_d,const unsigned int imgWidth,const unsigned int imgHeight)
{
	unsigned int tx =blockIdx.x*BLOCK_SIZE + threadIdx.x;
	unsigned int ty =blockIdx.y*BLOCK_SIZE + threadIdx.y;
	if(  tx < imgWidth &&  ty < imgHeight ){
		rhs_d[tx + imgWidth*ty]=extended_rhs_d[tx + 2*imgWidth*ty];
	}
}

__global__ void SpectralModulation(float2 *  extended_rhs_fft_d,const unsigned int extended_imgWidth,const unsigned int extended_imgHeight, const float finv_imgWidth,const float finv_imgHeight, const float dc)
{
	unsigned int tx =blockIdx.x*BLOCK_SIZE + threadIdx.x;
	unsigned int ty =blockIdx.y*BLOCK_SIZE + threadIdx.y;
	if(  tx < (extended_imgWidth/2 + 1)  &&  ty < extended_imgHeight ){
		if(tx + ty > 0)
		{
			float attenuation_factor = (finv_imgWidth*finv_imgHeight)/ (4.f -2.f*( cos(PI2*(float)tx*finv_imgWidth) + cos(PI2*((float)ty)*finv_imgHeight)));
			extended_rhs_fft_d[tx + ty*(extended_imgWidth/2 + 1) ].x *= attenuation_factor;
			extended_rhs_fft_d[tx + ty*(extended_imgWidth/2 + 1) ].y *= attenuation_factor;
		}
		else
		{
			extended_rhs_fft_d[0].x = 4.f*dc*(finv_imgWidth*finv_imgHeight);
			extended_rhs_fft_d[0].y = 0.f;
		}
	}
}


void GPUSolvers::FFTLaplaceSolver(float *  rhs_d, float * symmetric_extended_rhs_d, float2 * symmetric_extended_rhs_fft_d,  cufftHandle & fftPlanFwd, cufftHandle & fftPlanInv, const unsigned int array_width,const unsigned int array_height ,const float dc)
{
	unsigned int extended_array_width = 2*array_width;
	unsigned int extended_array_height = 2*array_height;

	//if(array_size!=array_width*array_height){
	//	array_size = array_width*array_height;
	//	
	//	if(symmetric_extended_rhs_d!=0){
	//		cudaFree(symmetric_extended_rhs_d);
	//	}
	//	cudaMalloc((void**)&symmetric_extended_rhs_d,extended_array_width*extended_array_height*sizeof(float));

	//	if(symmetric_extended_rhs_fft_d!=0){
	//		cudaFree(symmetric_extended_rhs_fft_d);
	//	}
	//	cudaMalloc((void**)&symmetric_extended_rhs_fft_d,(extended_array_height )*(extended_array_width/ 2 + 1) *sizeof(fcomplex));
	//}

	unsigned int  blocksW = (unsigned int) ceilf( (float) array_width / (float) BLOCK_SIZE );
    unsigned int  blocksH = (unsigned int) ceilf( (float) array_height /(float) BLOCK_SIZE );
    dim3 gridDim( blocksW, blocksH, 1 );
    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );

	//cudaEvent_t start;
	//cudaEventCreate(&start);
	//cudaEvent_t stop;
	//cudaEventCreate(&stop);
	//cudaEventRecord(start, NULL);

	ExtendedSymmetricRHS<<< gridDim, blockDim >>>( rhs_d, symmetric_extended_rhs_d, array_width, array_height);

	//cufftHandle fftPlanFwd, fftPlanInv;
 //   //printf("...creating R2C & C2R FFT plans for %i x %i\n", fftH, fftW);
 //   cufftSafeCall( cufftPlan2d(&fftPlanFwd, extended_array_height, extended_array_width, CUFFT_R2C) );
 //   cufftSafeCall( cufftPlan2d(&fftPlanInv, extended_array_height, extended_array_width, CUFFT_C2R) );

	cufftSafeCall( cufftExecR2C(fftPlanFwd, (cufftReal *)symmetric_extended_rhs_d, (cufftComplex *)symmetric_extended_rhs_fft_d) );
		
	unsigned int  ext_blocksW = (unsigned int) ceilf( (float) (extended_array_width / 2 + 1) / (float) BLOCK_SIZE );
    unsigned int  ext_blocksH = (unsigned int) ceilf( (float) extended_array_height /(float) BLOCK_SIZE );
    dim3 ext_gridDim( ext_blocksW, ext_blocksH, 1 );

	SpectralModulation<<< ext_gridDim, blockDim >>>(symmetric_extended_rhs_fft_d,extended_array_width,extended_array_height, 1.f/(float)extended_array_width,1.f/(float)extended_array_height,dc);

	cufftSafeCall( cufftExecC2R(fftPlanInv, (cufftComplex *)symmetric_extended_rhs_fft_d, (cufftReal *)symmetric_extended_rhs_d) );
	ExtractRHS<<< gridDim, blockDim >>>( rhs_d, symmetric_extended_rhs_d, array_width, array_height);

	//cudaEventRecord(stop, NULL);
	//cudaEventSynchronize(stop);
 //   float msecTotal = 0.0f;
 //   cudaEventElapsedTime(&msecTotal, start, stop);
	//printf("Time= %.5f msec \n",msecTotal);
}


//void GPUSolvers::FFTLaplaceSolver() // PASSED TEST
//{
//	int array_height = 359;
//	int array_width = 400;
//	const float dc = 57694.f;
//
//	unsigned int extended_array_width = 2*array_width;
//	unsigned int extended_array_height = 2*array_height;
//
//	if(symmetric_extended_rhs_d_size!=array_width*array_height){
//		if(rhs_d!=0){
//			cudaFree(rhs_d);
//		}
//		else{
//			cudaMalloc((void**)&rhs_d,array_width*array_height*sizeof(float));
//		}
//
//		if(symmetric_extended_rhs_d!=0){
//			cudaFree(symmetric_extended_rhs_d);
//		}
//		else{
//			cudaMalloc((void**)&symmetric_extended_rhs_d,extended_array_width*extended_array_height*sizeof(float));
//		}
//		if(symmetric_extended_rhs_h!=0){
//			delete symmetric_extended_rhs_h;
//		}
//		else{
//			symmetric_extended_rhs_h = new float[extended_array_width*extended_array_height*sizeof(float)];
//		}
//		if(symmetric_extended_rhs_fft_d!=0){
//			cudaFree(symmetric_extended_rhs_fft_d);
//		}
//		else{
//			cudaMalloc((void**)&symmetric_extended_rhs_fft_d,(extended_array_height )*(extended_array_width/ 2 + 1) *sizeof(fcomplex));
//		}
//		if(symmetric_extended_rhs_fft_h!=0){
//			delete symmetric_extended_rhs_fft_h;
//		}
//		else{
//			symmetric_extended_rhs_fft_h = new fcomplex[(extended_array_height)*(extended_array_width/ 2 + 1) *sizeof(fcomplex)];
//		}
//		if(rhs_h!=0){
//			delete rhs_h;
//		}
//		else{
//			rhs_h = new float[array_width*array_height*sizeof(float)];
//		}
//	}
//
//	FileIO::readInputdf("input.txt",rhs_h,array_height*array_width);
//
//	printf("input [0] = %g \n",rhs_h[0]);
//	printf("input [1] = %g \n",rhs_h[1]);
//	printf("input [width - 1] = %g \n",rhs_h[array_width-1]);
//	printf("input [width] = %g \n",rhs_h[array_width]);
//
//	cudaMemcpy(rhs_d,rhs_h,array_width*array_height*sizeof(float),cudaMemcpyHostToDevice);
//
//	unsigned int  blocksW = (unsigned int) ceilf( (float) array_width / (float) BLOCK_SIZE );
//    unsigned int  blocksH = (unsigned int) ceilf( (float) array_height /(float) BLOCK_SIZE );
//    dim3 gridDim( blocksW, blocksH, 1 );
//    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );
//
//	ExtendedSymmetricRHS<<< gridDim, blockDim >>>( rhs_d, symmetric_extended_rhs_d, array_width, array_height);
//
//	cudaMemcpy(symmetric_extended_rhs_h,symmetric_extended_rhs_d,extended_array_width*extended_array_height*sizeof(float),cudaMemcpyDeviceToHost);
//
//	printf("extended [0] = %g \n",symmetric_extended_rhs_h[0]);
//	printf("extended [1] = %g \n",symmetric_extended_rhs_h[1]);
//	printf("extended [array_width-2] = %g \n",symmetric_extended_rhs_h[array_width-2]);
//	printf("extended [array_width-1] = %g \n",symmetric_extended_rhs_h[array_width-1]);
//	printf("extended [array_width] = %g \n",symmetric_extended_rhs_h[array_width]);
//	printf("extended [array_width+1] = %g \n",symmetric_extended_rhs_h[array_width+1]);
//	printf("extended [array_width+2] = %g \n",symmetric_extended_rhs_h[array_width+2]);
//	printf("extended [extended_array_width-1] = %g \n",symmetric_extended_rhs_h[extended_array_width-1]);
//	printf("extended [extended_array_width] = %g \n",symmetric_extended_rhs_h[extended_array_width]);
//
//	cufftHandle fftPlanFwd, fftPlanInv;
//    //printf("...creating R2C & C2R FFT plans for %i x %i\n", fftH, fftW);
//    cufftSafeCall( cufftPlan2d(&fftPlanFwd, extended_array_height, extended_array_width, CUFFT_R2C) );
//    cufftSafeCall( cufftPlan2d(&fftPlanInv, extended_array_height, extended_array_width, CUFFT_C2R) );
//	cufftSafeCall( cufftExecR2C(fftPlanFwd, (cufftReal *)symmetric_extended_rhs_d, (cufftComplex *)symmetric_extended_rhs_fft_d) );
//		
//	cudaMemcpy(symmetric_extended_rhs_fft_h,symmetric_extended_rhs_fft_d,(extended_array_height)*(extended_array_width / 2 + 1) *sizeof(fcomplex),cudaMemcpyDeviceToHost);
//
//	printf("extended fft [0] = (%g,%g) \n",symmetric_extended_rhs_fft_h[0].x,symmetric_extended_rhs_fft_h[0].y);
//	printf("extended fft [1] = (%g,%g) \n",symmetric_extended_rhs_fft_h[1].x,symmetric_extended_rhs_fft_h[1].y);
//	printf("extended fft [array_width-2] = (%g,%g) \n",symmetric_extended_rhs_fft_h[array_width-2].x,symmetric_extended_rhs_fft_h[array_width-2].y);
//	printf("extended fft [array_width-1] = (%g,%g) \n",symmetric_extended_rhs_fft_h[array_width-1].x,symmetric_extended_rhs_fft_h[array_width-1].y);
//	printf("extended fft [array_width] = (%g,%g) \n",symmetric_extended_rhs_fft_h[array_width].x,symmetric_extended_rhs_fft_h[array_width].y);
//	printf("extended fft [array_width+1] = (%g,%g) \n",symmetric_extended_rhs_fft_h[array_width+1].x,symmetric_extended_rhs_fft_h[array_width+1].y);
//	printf("extended fft [array_width+2] = (%g,%g) \n",symmetric_extended_rhs_fft_h[array_width+2].x,symmetric_extended_rhs_fft_h[array_width+2].y);
//	printf("extended fft [extended_array_width-1] = (%g,%g) \n",symmetric_extended_rhs_fft_h[extended_array_width-1].x,symmetric_extended_rhs_fft_h[extended_array_width-1].y);
//	printf("extended fft [extended_array_width] = (%g,%g) \n",symmetric_extended_rhs_fft_h[extended_array_width].x,symmetric_extended_rhs_fft_h[extended_array_width].y);
//
//
//	unsigned int  ext_blocksW = (unsigned int) ceilf( (float) (extended_array_width / 2 + 1) / (float) BLOCK_SIZE );
//    unsigned int  ext_blocksH = (unsigned int) ceilf( (float) extended_array_height /(float) BLOCK_SIZE );
//    dim3 ext_gridDim( ext_blocksW, ext_blocksH, 1 );
//
//	//cudaEvent_t start;
//	//cudaEventCreate(&start);
//	//cudaEvent_t stop;
//	//cudaEventCreate(&stop);
//	//cudaEventRecord(start, NULL);
//	
//	SpectralModulation<<< ext_gridDim, blockDim >>>(symmetric_extended_rhs_fft_d,extended_array_width,extended_array_height, 1.f/(float)extended_array_width,1.f/(float)extended_array_height,dc);
//	
//	cudaMemcpy(symmetric_extended_rhs_fft_h,symmetric_extended_rhs_fft_d,(extended_array_height)*(extended_array_width / 2 + 1) *sizeof(fcomplex),cudaMemcpyDeviceToHost);
//
//	printf("extended modulated fft [0] = (%g,%g) \n",symmetric_extended_rhs_fft_h[0].x,symmetric_extended_rhs_fft_h[0].y);
//	printf("extended modulated fft [1] = (%g,%g) \n",symmetric_extended_rhs_fft_h[1].x,symmetric_extended_rhs_fft_h[1].y);
//	printf("extended modulated fft [array_width-2] = (%g,%g) \n",symmetric_extended_rhs_fft_h[array_width-2].x,symmetric_extended_rhs_fft_h[array_width-2].y);
//	printf("extended modulated fft [array_width-1] = (%g,%g) \n",symmetric_extended_rhs_fft_h[array_width-1].x,symmetric_extended_rhs_fft_h[array_width-1].y);
//	printf("extended modulated fft [array_width] = (%g,%g) \n",symmetric_extended_rhs_fft_h[array_width].x,symmetric_extended_rhs_fft_h[array_width].y);
//	printf("extended modulated fft [array_width+1] = (%g,%g) \n",symmetric_extended_rhs_fft_h[array_width+1].x,symmetric_extended_rhs_fft_h[array_width+1].y);
//	printf("extended modulated fft [array_width+2] = (%g,%g) \n",symmetric_extended_rhs_fft_h[array_width+2].x,symmetric_extended_rhs_fft_h[array_width+2].y);
//	printf("extended modulated fft [extended_array_width-1] = (%g,%g) \n",symmetric_extended_rhs_fft_h[extended_array_width-1].x,symmetric_extended_rhs_fft_h[extended_array_width-1].y);
//	printf("extended modulated fft [extended_array_width] = (%g,%g) \n",symmetric_extended_rhs_fft_h[extended_array_width].x,symmetric_extended_rhs_fft_h[extended_array_width].y);
//
//
//	//cudaEventRecord(stop, NULL);
//	//cudaEventSynchronize(stop);
// //   float msecTotal = 0.0f;
// //   cudaEventElapsedTime(&msecTotal, start, stop);
//	//printf("Time= %.5f msec \n",msecTotal);
//
//
//	cufftSafeCall( cufftExecC2R(fftPlanInv, (cufftComplex *)symmetric_extended_rhs_fft_d, (cufftReal *)symmetric_extended_rhs_d) );
//	cudaMemcpy(symmetric_extended_rhs_h,symmetric_extended_rhs_d,extended_array_width*extended_array_height*sizeof(float),cudaMemcpyDeviceToHost);
//	
//	printf("extended [0] = %g \n",symmetric_extended_rhs_h[0]);
//	printf("extended [1] = %g \n",symmetric_extended_rhs_h[1]);
//	printf("extended [array_width-2] = %g \n",symmetric_extended_rhs_h[array_width-2]);
//	printf("extended [array_width-1] = %g \n",symmetric_extended_rhs_h[array_width-1]);
//	printf("extended [array_width] = %g \n",symmetric_extended_rhs_h[array_width]);
//	printf("extended [array_width+1] = %g \n",symmetric_extended_rhs_h[array_width+1]);
//	printf("extended [array_width+2] = %g \n",symmetric_extended_rhs_h[array_width+2]);
//	printf("extended [extended_array_width-1] = %g \n",symmetric_extended_rhs_h[extended_array_width-1]);
//	printf("extended [extended_array_width] = %g \n",symmetric_extended_rhs_h[extended_array_width]);
//
//	ExtractRHS<<< gridDim, blockDim >>>( rhs_d, symmetric_extended_rhs_d, array_width, array_height);
//	cudaMemcpy(rhs_h,rhs_d,array_width*array_height*sizeof(float),cudaMemcpyDeviceToHost);
//
//	printf("res [0] = %g \n",rhs_h[0]);
//	printf("res [1] = %g \n",rhs_h[1]);
//	printf("res [array_width-2] = %g \n",rhs_h[array_width-2]);
//	printf("res [array_width-1] = %g \n",rhs_h[array_width-1]);
//	printf("res [array_width] = %g \n",rhs_h[array_width]);
//	printf("res [array_width+1] = %g \n",rhs_h[array_width+1]);
//	printf("res [array_width+2] = %g \n",rhs_h[array_width+2]);
//
//	FileIO::writeOutputff("output.txt",rhs_h,array_height*array_width);
//}

//void GPUSolvers::FFTLaplaceSolver(float *  rhs_d, const unsigned int array_width,const unsigned int array_height)
//{
//	unsigned int extended_array_width = 2*array_width;
//	unsigned int extended_array_height = 2*array_height;
//
//	if(symmetric_extended_rhs_d_size!=array_width*array_height){
//		if(symmetric_extended_rhs_d!=0){
//			cudaFree(symmetric_extended_rhs_d);
//		}
//		else{
//			cudaMalloc((void**)&symmetric_extended_rhs_d,extended_array_width*extended_array_height*sizeof(float));
//		}
//		if(symmetric_extended_rhs_fft_d!=0){
//			cudaFree(symmetric_extended_rhs_fft_d);
//		}
//		else{
//			cudaMalloc((void**)&symmetric_extended_rhs_fft_d,extended_array_height*(extended_array_width / 2 + 1) *sizeof(fcomplex));
//		}
//	}
//
//
//
//	cudaMemcpy(rhs_h,rhs_d,array_width*array_height*sizeof(float),cudaMemcpyDeviceToHost);
//
//	printf("Entry 0 = %f", rhs_h[0]);
//
//	unsigned int  blocksW = (unsigned int) ceilf( (float) array_width / (float) BLOCK_SIZE );
//    unsigned int  blocksH = (unsigned int) ceilf( (float) array_height /(float) BLOCK_SIZE );
//    dim3 gridDim( blocksW, blocksH, 1 );
//    dim3 blockDim( BLOCK_SIZE, BLOCK_SIZE, 1 );
//
//	ExtendedSymmetricRHS<<< gridDim, blockDim >>>( rhs_d, symmetric_extended_rhs_d, array_width, array_height);
//
//	cufftHandle fftPlanFwd, fftPlanInv;
//    //printf("...creating R2C & C2R FFT plans for %i x %i\n", fftH, fftW);
//    cufftSafeCall( cufftPlan2d(&fftPlanFwd, extended_array_height, extended_array_width, CUFFT_R2C) );
//    cufftSafeCall( cufftPlan2d(&fftPlanInv, extended_array_height, extended_array_width, CUFFT_C2R) );
//	cufftSafeCall( cufftExecR2C(fftPlanFwd, (cufftReal *)symmetric_extended_rhs_d, (cufftComplex *)symmetric_extended_rhs_fft_d) );
//
//	unsigned int  ext_blocksW = (unsigned int) ceilf( (float) (extended_array_width / 2 + 1) / (float) BLOCK_SIZE );
//    unsigned int  ext_blocksH = (unsigned int) ceilf( (float) extended_array_height /(float) BLOCK_SIZE );
//    dim3 ext_gridDim( ext_blocksW, ext_blocksH, 1 );
//
//	//cudaEvent_t start;
//	//cudaEventCreate(&start);
//	//cudaEvent_t stop;
//	//cudaEventCreate(&stop);
//	//cudaEventRecord(start, NULL);
//	
//	SpectralInversion<<< ext_gridDim, blockDim >>>(symmetric_extended_rhs_fft_d,extended_array_width,extended_array_height, 1.f/(float)extended_array_width,1.f/(float)extended_array_height);
//	
//	//cudaEventRecord(stop, NULL);
//	//cudaEventSynchronize(stop);
// //   float msecTotal = 0.0f;
// //   cudaEventElapsedTime(&msecTotal, start, stop);
//	//printf("Time= %.5f msec \n",msecTotal);
//
//
//	cufftSafeCall( cufftExecC2R(fftPlanInv, (cufftComplex *)symmetric_extended_rhs_fft_d, (cufftReal *)symmetric_extended_rhs_d) );
//
//	cudaMemcpy(symmetric_extended_rhs_fft_h,symmetric_extended_rhs_fft_d,extended_array_height*(extended_array_width / 2 + 1) *sizeof(fcomplex),cudaMemcpyDeviceToHost);
//
//	ExtractRHS<<< gridDim, blockDim >>>( rhs_d, symmetric_extended_rhs_d, array_width, array_height);
//}

