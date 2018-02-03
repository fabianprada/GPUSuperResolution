#ifndef GL_TOOLS_INCLUDED
#define GL_TOOLS_INCLUDED

#include "GL\glew.h"
#include "GL\glut.h"
#include <stdio.h>

#define INCLUDE_CUDA 1
#if INCLUDE_CUDA
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
#endif

namespace GLTools
{

	void CreateFramebuffer( GLuint& framebuffer, GLuint colorAttachment0, GLuint depthAttachment );
	void DeleteFramebuffer( GLuint& framebuffer );

	void CreateVBO( GLuint& bufferID, unsigned int size );
	void DeleteVBO( GLuint& bufferID );

	void CreatePBO( GLuint& bufferID, size_t size );
	void DeletePBO( GLuint& bufferID );

	void CreateTexture( GLuint& texture,const unsigned int width,const unsigned int height, unsigned char * u_values = NULL ,float * f_values = NULL,const unsigned int input_num_channels = 4, const unsigned int texture_num_channels = 4,GLenum texture_type = GL_UNSIGNED_BYTE,GLenum texture_boundary = GL_MIRRORED_REPEAT,GLenum texture_filter =GL_NEAREST);
	void DeleteTexture( GLuint& texture );

	void CreateDepthBuffer( GLuint& depthBuffer, unsigned int width, unsigned int height );
	void DeleteDepthBuffer( GLuint& depthBuffer );

	void TransformFloat4ArrayToGLColorBuffer_CUDA(cudaGraphicsResource_t& cuda_resource, float4 * float4_buffer, unsigned char * uchar4_buffer, const unsigned int width,const unsigned int height);
	void TransformFloat2ArrayToGLVectorBuffer_CUDA(cudaGraphicsResource_t& cuda_resource, float2 * vector_buffer, const unsigned int imgWidth,const unsigned int imgHeigh, bool normalize_vector_visualization,bool clamp_for_visualization,const float clamping_threshold,const float field_amplification);
	void Transform2FloatArrayToGLVectorBuffer_CUDA(cudaGraphicsResource_t& cuda_resource, float * vector0_buffer,float * vector1_buffer, const unsigned int imgWidth,const unsigned int imgHeigh, bool normalize_vector_visualization,bool clamp_for_visualization,const float clamping_threshold,const float field_amplification);
	void NormalizeVectorFieldByMaxima(float2 * vector_buffer, float * norm_buffer, int width, int height);
	void GeneralVectorScaling(float2 * vector_buffer, float * norm_buffer, int width, int height);
	// Links a OpenGL texture object to a CUDA resource that can be used in the CUDA kernel.

#if INCLUDE_CUDA
	void CreateCUDAImageResource( cudaGraphicsResource_t& cudaResource, GLuint GLtexture, cudaGraphicsMapFlags mapFlags );
	void CreateCUDABufferResource( cudaGraphicsResource_t& cudaResource, GLuint GLbuffer, cudaGraphicsMapFlags mapFlags );
	void DeleteCUDAResource( cudaGraphicsResource_t& cudaResource );
#endif
}

#endif// GL_TOOLS_INCLUDED