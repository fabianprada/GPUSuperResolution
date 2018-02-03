#include <cufft.h>
namespace GPUSolvers{
  void FFTLaplaceSolver(float *  rhs_d,  float * symmetric_extended_rhs_d, float2 * symmetric_extended_rhs_fft_d, cufftHandle & fftPlanFwd, cufftHandle & fftPlanInv, const unsigned int array_width,const unsigned int array_height,const float dc); 
}