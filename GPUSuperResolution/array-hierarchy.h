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
#include <cufft.h>

class Array_Hierarchy : public CUDA_Advection_Object {
public:
	Array_Hierarchy(int p_hierarchy_levels , int p_height, int p_width, unsigned char * p_color, Canvas_Rectangle * p_canvas, Canvas_Rectangle * p_reference_canvas) : CUDA_Advection_Object(p_height,p_width,p_color,p_canvas){
		hierarchy_levels = p_hierarchy_levels;
		active_level = 2;
		padding_factor_threshold = 0.f;

		level_absolute_coordinate = new Point<2>[hierarchy_levels];
		level_absolute_scale = new double[hierarchy_levels];
		updated_level = new bool[hierarchy_levels];

		double current_scale =1.f;
		for(int i=0; i<hierarchy_levels; i++){
			updated_level[i] =false;
			level_absolute_coordinate[i]=Point<2>();
			level_absolute_scale[i]= current_scale;
			current_scale/=2.f;
		}

		reference_canvas = p_reference_canvas;

		//advection_steps =10;
	}

	int hierarchy_levels;
	float padding_factor_threshold;
	int active_level;

	Point<2> * level_absolute_coordinate;
	double * level_absolute_scale;
	bool * updated_level;

	void UpdateHierachy();
	void UpdateLevel(const int num_iter, const int level_number);
	void SamplePreviousLevel_CUDA(const int level_number);
	void AdvectCurrentLevel_CUDA(const int level_number);
	//void SampleAdvectionField_CUDA();
	//void BackwardAdvection_CUDA(const int num_iter, const int level_number);
	void SetTexture_CUDA(cudaGraphicsResource_t& p_texture_resource, unsigned char * uchar4_buffer_d, float4 * flaot4_buffer_d);
	void SetTextureToLevel(const int level_number,cudaGraphicsResource_t& p_texture_resource);
	
	cudaArray ** hierarchical_color_field_array;
	Canvas_Rectangle * reference_canvas;
	void setupOpenGL();
	void drawOpenGL(Scene * scene = 0);
	void HandleKeyboardEvent(Prompt  & prompt, Scene * scene = 0);
	void drawROI(const Point<2> p_roi_absolute_coordinate, const double p_roi_absolute_scale);
};
