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


#include "array.h"

template<int dimension>
void Array2D<dimension>::setupOpenGL()
{
	setupTexture(pixel_texture_buffer);
}

template<int dimension>
void Array2D<dimension>::setupTexture(GLuint & texture_buffer)
{
	if (!glIsTexture(texture_buffer))
	{
		glGenTextures(1, &texture_buffer);
		glBindTexture(GL_TEXTURE_2D, texture_buffer);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	}
	
}

template<int dimension>
double Array2D<dimension>::TransformValueForVisualization(double x)
{
	if (abs_visualization)
	{
		x = abs(x);
	}
	if (active_percentile_visualization > -1)
	{
		double scale_factor = std::max< double >(x, statistics->value_linf_norm_percentiles[active_percentile_visualization]);
		x /= scale_factor;
	}
	return x;
	
}
template<int dimension>
void Array2D<dimension>::setTextureValues(GLuint & texture_buffer)
{
	glBindTexture(GL_TEXTURE_2D, texture_buffer);

	unsigned char * img_data = new unsigned char[height*width*dimension];
	int position = 0; /// parallel for !!!

#pragma omp parallel for
	for (int j = 0; j < height; j++){
		int write_position = j*width*dimension;
		int read_position =  j*width;
		for (int i = 0; i < width; i++){
			Point<dimension> color_value = entries[read_position];
			for (int k = 0; k < dimension; k++){
				int color = static_cast<int>(round(color_value[k] * 255.f));
				unsigned char ucolor;
				if (color < 0)
					ucolor = 0;
				else if (color>255)
					ucolor = 255;
				else
					ucolor = static_cast<unsigned char>(color);
				img_data[write_position] = ucolor;
				write_position++;
			}
			read_position++;
		}
	}
	if (dimension == 1){
		glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, width, height, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, img_data);
	}
	else if (dimension == 3){
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, img_data);
	}
	else if (dimension == 4){
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	}
	else{
		printf("Unkown Texture Conversion");
	}
	delete[] img_data;

	pixel_texture_buffer_updated = true;
}

template<int dimension>
void Array2D<dimension>::setTextureValues(GLuint & texture_buffer, const int corner_w, const int corner_h, const int region_w, const int region_h)
{
	const int last_pos_h = corner_h + region_h;
	const int last_pos_w = corner_w + region_w;
	if (last_pos_w > width || last_pos_h > height){
		printf("Unvalid Dimensions. Unable Set Texture Values !! \n");
	}

	glBindTexture(GL_TEXTURE_2D, texture_buffer);
	unsigned char * img_data = new unsigned char[region_w*region_h*dimension];

//#pragma omp parallel for

	for (int j = 0; j < region_h; j++){
		int write_position = j*region_w*dimension;
		int read_position = (j + corner_h)*width + corner_w;
		for (int i = 0; i < region_w; i++){
			Point<dimension> color_value = entries[read_position];
			for (int k = 0; k < dimension; k++){
				int color = static_cast<int>(round(color_value[k] * 255.f));
				unsigned char ucolor;
				if (color < 0)
					ucolor = 0;
				else if (color>255)
					ucolor = 255;
				else
					ucolor = static_cast<unsigned char>(color);
				img_data[write_position] = ucolor;
				write_position++;
			}
			read_position++;
		}
	}

	if (dimension == 1){
		glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, region_w, region_h, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, img_data);
	}
	else if (dimension == 3){
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, region_w, region_h, 0, GL_RGB, GL_UNSIGNED_BYTE, img_data);
	}
	else if (dimension == 4){
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, region_w, region_h, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	}
	else{
		printf("Unkown Texture Conversion");
	}
	delete[] img_data;
}

template<int dimension>
void Array2D<dimension>::drawOpenGL(Scene * scene = 0)
{
	if (!pixel_texture_buffer_updated){
		setTextureValues(pixel_texture_buffer);
	}
	glEnable(GL_TEXTURE_2D);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glBindTexture(GL_TEXTURE_2D, pixel_texture_buffer);

	glBegin(GL_QUADS);
	Point<3> vertex_position = canvas->reference_vertex;
	glTexCoord2d(0.f, 0.f);
	glVertex3d(vertex_position[0], vertex_position[1], vertex_position[2]);

	vertex_position += canvas->axis[0] * canvas->edges_lenght[0];
	glTexCoord2d(1.f, 0.f);
	glVertex3d(vertex_position[0], vertex_position[1], vertex_position[2]);

	vertex_position += canvas->axis[1] * canvas->edges_lenght[1];
	glTexCoord2d(1.f, 1.f);
	glVertex3d(vertex_position[0], vertex_position[1], vertex_position[2]);

	vertex_position -= canvas->axis[0] * canvas->edges_lenght[0];
	glTexCoord2d(0.f, 1.f);
	glVertex3d(vertex_position[0], vertex_position[1], vertex_position[2]);

	glEnd();

	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisable(GL_TEXTURE_2D);
}


