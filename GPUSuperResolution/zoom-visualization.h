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


#ifndef ZOOM_VISUALIZATION_INCLUDED
#define ZOOM_VISUALIZATION_INCLUDED

#include "graphic-object.h"

class ZoomInObject : public GraphicObject{
public:
	virtual void drawROI(const Point<2> roi_absolute_coordinate, const double roi_absolute_scale) = 0;
};

enum SamplingDomain { GLOBAL_DOMAIN, LOCAL_DOMAIN, HYBRID_DOMAIN};
enum VectorFieldMode { MEAN_GRADIENT_FIELD_LINEAR,MEAN_GRADIENT_FIELD, ADVECTION_VECTOR_FIELD, NONE_VECTOR_FIELD};
enum GradientAdvectionMode { PARENT_STATE_SAMPLED_GRADIENT, PARENT_STATE_ACCUMULATED_GRADIENT};
enum NormalizationMode { UNIFORM_NORMALIZATION, MAXIMA_NORMALIZATION, NONE_NORMALIZATION };

class ZoomVisualization : public GraphicObject{
public:
	ZoomVisualization(int p_height, int p_width, unsigned char * p_color, ZoomInObject * p_zoom_object, Canvas_Rectangle * p_reference_canvas,Point<2> p_roi_absolute_coordinate = Point<2>(0.125f,0.125f), double p_roi_absolute_scale = 0.499f, double p_roi_padding_factor = 0.0625f){//int p_num_levels,

		sprintf(object_name, "ZOOM VISUALIZATION");
		is_active = true;

		height = p_height;
		width = p_width;
		size = height*width;

		color = p_color;

		roi_absolute_coordinate = p_roi_absolute_coordinate;
		roi_absolute_scale = p_roi_absolute_scale;
		roi_padding_factor = p_roi_padding_factor;

		zoom_object = p_zoom_object;

		reference_canvas = p_reference_canvas;
	}

	int height;
	int width;
	int size; 

	unsigned char * color;
	Canvas_Rectangle * reference_canvas;

	Point<2> roi_absolute_coordinate;
	double roi_absolute_scale;
	double roi_padding_factor;

	GLuint color_texture_buffer;

	ZoomInObject * zoom_object;


	void setupOpenGL();
	void drawOpenGL(Scene * scene = 0);
	void drawBox(Canvas_Rectangle * p_canvas, Point<2> p_canvas_coordinates, double p_canvas_scale, Point<3> p_color = Point<3>(), double p_width = 0.01f);
	void HandleMouseEvent(MouseEvent  & mouse_event, Scene * scene = 0);
	void HandleKeyboardEvent(Prompt  & prompt, Scene * scene = 0);
	
};

#endif // ZOOM_VISUALIZATION_INCLUDED