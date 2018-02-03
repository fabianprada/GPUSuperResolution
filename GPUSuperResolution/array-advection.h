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

class Array_Advection : public CUDA_Advection_Object{
public:
	Array_Advection(int p_height, int p_width, unsigned char * p_color, Canvas_Rectangle * p_canvas, Point<2> p_absolute_coordinate = Point<2>(), double p_absolute_scale = 1.f) : CUDA_Advection_Object(p_height,p_width,p_color,p_canvas){
	
		roi_absolute_coordinate = p_absolute_coordinate;
		roi_absolute_scale = p_absolute_scale;
	}

	void setupOpenGL();
	void drawOpenGL(Scene * scene = 0);
	void HandleKeyboardEvent(Prompt  & prompt, Scene * scene = 0);
	void drawROI(const Point<2> p_roi_absolute_coordinate, const double p_roi_absolute_scale);
};
