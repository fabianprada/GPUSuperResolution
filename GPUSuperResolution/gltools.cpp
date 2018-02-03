#include "gltools.h"
#include <math.h>

// Create a pixel buffer object
void GLTools::CreateVBO( GLuint& bufferID, unsigned int size)
{
    // Make sure the buffer doesn't already exist
    DeleteVBO( bufferID );

    glGenBuffers( 1, &bufferID );
    glBindBuffer( GL_ARRAY_BUFFER, bufferID );
    glBufferData( GL_ARRAY_BUFFER, size, NULL, GL_STREAM_DRAW );

    glBindBuffer( GL_ARRAY_BUFFER, 0 );
}

void GLTools::DeleteVBO(  GLuint& bufferID )
{
    if ( bufferID != 0 )
    {
        glDeleteBuffers( 1, &bufferID );
        bufferID = 0;
    }
}

// Create a pixel buffer object
void GLTools::CreatePBO( GLuint& bufferID, size_t size )
{
    // Make sure the buffer doesn't already exist
    DeletePBO( bufferID );

    glGenBuffers( 1, &bufferID );
    glBindBuffer( GL_PIXEL_UNPACK_BUFFER, bufferID );
    glBufferData( GL_PIXEL_UNPACK_BUFFER, size, NULL, GL_STREAM_DRAW );

    glBindBuffer( GL_PIXEL_UNPACK_BUFFER, 0 );
}

void GLTools::DeletePBO(  GLuint& bufferID )
{
    if ( bufferID != 0 )
    {
        glDeleteBuffers( 1, &bufferID );
        bufferID = 0;
    }
}

unsigned char NormalFloatToUnsigned( float x)
{
	int i_color = static_cast<int>(floor(x*255.f));
	if (i_color < 0)
		i_color = 0;
	else if (i_color>255)
		i_color = 255;
	return static_cast<unsigned char>(i_color);
}



// Create a texture resource for rendering to.
void GLTools::CreateTexture( GLuint& texture,const unsigned int width,const unsigned int height, unsigned char * u_values, float * f_values, const unsigned int input_num_channels, const unsigned int texture_num_channels, GLenum texture_type,GLenum texture_boundary,GLenum texture_filter)
{
    // Make sure we don't already have a texture defined here
    DeleteTexture( texture );

    glGenTextures( 1, &texture );
    glBindTexture( GL_TEXTURE_2D, texture );

    // set basic parameters
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, texture_boundary);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, texture_boundary);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, texture_filter);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, texture_filter);

    // Create texture data (4-component unsigned byte)

	if(u_values != NULL && f_values != NULL){
		printf("UNEXPECTED CONDITION!!\n");
	}

	if(texture_type == GL_FLOAT)
	{
		float * buffer_values;
		if(u_values!= NULL){
			buffer_values = new float[height*width*texture_num_channels];
			for(int i=0; i<height*width; i++){
				if(texture_num_channels == input_num_channels){
					for(int k =0; k < input_num_channels; k++){
						buffer_values[i*input_num_channels + k] = ((float)u_values[i*input_num_channels + k])/255.f;
					}
				}
				else if(texture_num_channels == 1 &&  input_num_channels == 3){
					buffer_values[i] = (((float)u_values[i*3])*0.3f + ((float)u_values[i*3 + 1])*0.59f + ((float)u_values[i*3 + 2])*0.11f) /255.f;
				}
				else if(texture_num_channels == 1 &&  input_num_channels == 4){
					buffer_values[i] = (((float)u_values[i*4])*0.3f + ((float)u_values[i*4 + 1])*0.59f + ((float)u_values[i*4 + 2])*0.11f) /255.f;
				}
				else{
					printf("Undefined Dimension Transformation!!\n");
				}
			}
		}
		else{
			buffer_values =f_values;
		}

		if(texture_num_channels ==4){
			glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA,texture_type,buffer_values);
		}
		else if(texture_num_channels ==3){
			glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB,texture_type,buffer_values);
		}
		else if(texture_num_channels ==2){
			glTexImage2D( GL_TEXTURE_2D, 0, GL_RG, width, height, 0, GL_RG,texture_type,buffer_values);
		}
		else if(texture_num_channels ==1){
			glTexImage2D( GL_TEXTURE_2D, 0, GL_R8, width, height, 0, GL_LUMINANCE,texture_type,buffer_values);
		}
		else{
			printf("UNEXPECTED TEXTURE DEPTH\n");
		}

		if(buffer_values != f_values){
			delete buffer_values;
		}
	}

	if(texture_type == GL_UNSIGNED_BYTE)
	{
		unsigned char * buffer_values;
		if(f_values!= NULL){
			buffer_values = new unsigned char[height*width*texture_num_channels];
			for(int i=0; i<height*width;i++){

				if(texture_num_channels == input_num_channels){
					for(int k =0; k < input_num_channels; k++){
					buffer_values[i*input_num_channels + k] = NormalFloatToUnsigned(f_values[i*input_num_channels + k]);
					}
				}
				else if(texture_num_channels == 1 &&  input_num_channels == 3){
					buffer_values[i] = NormalFloatToUnsigned((f_values[i*3])*0.3f + (f_values[i*3 + 1])*0.59f + (f_values[i*3 + 2])*0.11f);
				}
				else if(texture_num_channels == 1 &&  input_num_channels == 4){
					buffer_values[i] = NormalFloatToUnsigned((f_values[i*4])*0.3f + (f_values[i*4 + 1])*0.59f + (f_values[i*4 + 2])*0.11f);
				}
				else{
					printf("Undefined Dimension Transformation!!\n");
				}
			}
		}
		else{
			buffer_values = u_values;
		}

		if(texture_num_channels ==4){
			glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA,texture_type,buffer_values);
		}
		else if(texture_num_channels ==3){
			glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB,texture_type,buffer_values);
		}
		else if(texture_num_channels ==2){
			glTexImage2D( GL_TEXTURE_2D, 0, GL_RG, width, height, 0, GL_RG,texture_type,buffer_values);
		}
		else if(texture_num_channels ==1){
			glTexImage2D( GL_TEXTURE_2D, 0, GL_LUMINANCE, width, height, 0, GL_LUMINANCE,texture_type,buffer_values);
		}
		else{
			printf("UNEXPECTED TEXTURE DEPTH\n");
		}

		if(buffer_values != u_values){
			delete buffer_values;
		}
	}
	
    // Unbind the texture
    glBindTexture( GL_TEXTURE_2D, 0 );
}

void GLTools::DeleteTexture( GLuint& texture )
{
    if ( texture != 0 )
    {
        glDeleteTextures(1, &texture );
        texture = 0;
    }
}

void GLTools::CreateDepthBuffer( GLuint& depthBuffer, unsigned int width, unsigned int height )
{
    // Delete the existing depth buffer if there is one.
    DeleteDepthBuffer( depthBuffer );

    glGenRenderbuffers( 1, &depthBuffer );
    glBindRenderbuffer( GL_RENDERBUFFER, depthBuffer );

    glRenderbufferStorage( GL_RENDERBUFFER, GL_DEPTH_COMPONENT, width, height );

    // Unbind the depth buffer
    glBindRenderbuffer( GL_RENDERBUFFER, 0 );
}

void GLTools::DeleteDepthBuffer( GLuint& depthBuffer )
{
    if ( depthBuffer != 0 )
    {
        glDeleteRenderbuffers( 1, &depthBuffer );
        depthBuffer = 0;
    }
}

void GLTools::CreateFramebuffer( GLuint& framebuffer, GLuint colorAttachment0, GLuint depthAttachment )
{
    // Delete the existing framebuffer if it exists.
    DeleteFramebuffer( framebuffer );

    glGenFramebuffers( 1, &framebuffer );
    glBindFramebuffer( GL_FRAMEBUFFER, framebuffer );

    glFramebufferTexture2D( GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, colorAttachment0, 0 );
    glFramebufferRenderbuffer( GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthAttachment );

    // Check to see if the frame buffer is valid
    GLenum fboStatus = glCheckFramebufferStatus( GL_FRAMEBUFFER );
    if ( fboStatus != GL_FRAMEBUFFER_COMPLETE )
    {
        printf("ERROR: Incomplete framebuffer status.\n");
    }

    // Unbind the frame buffer
    glBindFramebuffer( GL_FRAMEBUFFER, 0 );
}

void GLTools::DeleteFramebuffer( GLuint& framebuffer )
{
    if ( framebuffer != 0 )
    {
        glDeleteFramebuffers( 1, &framebuffer );
        framebuffer = 0;
    }
}

void GLTools::CreateCUDAImageResource( cudaGraphicsResource_t& cudaResource, GLuint GLtexture, cudaGraphicsMapFlags mapFlags )
{
    // Map the GL texture resource with the CUDA resource
    cudaGraphicsGLRegisterImage( &cudaResource, GLtexture, GL_TEXTURE_2D, mapFlags );
}

void GLTools::CreateCUDABufferResource( cudaGraphicsResource_t& cudaResource, GLuint GLbuffer, cudaGraphicsMapFlags mapFlags )
{
    // Map the GL texture resource with the CUDA resource
    cudaGraphicsGLRegisterBuffer( &cudaResource, GLbuffer, mapFlags );
}

void GLTools::DeleteCUDAResource( cudaGraphicsResource_t& cudaResource )
{
    if ( cudaResource != 0 )
    {
        cudaGraphicsUnregisterResource( cudaResource );
        cudaResource = 0;
    }
}