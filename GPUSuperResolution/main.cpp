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


#include <cstdlib>
#include <iostream>
#include <string>
#include "window.h"
#include "image.h"
#include "zoom-visualization.h"
#include "array-advection.h"
#include "array-hierarchy.h"
#include "filter.h"
#include "extension.h"
// CUDA headers
#include "cpu-solver.h"

#define PI 3.1415926
#define WINDOW_HEIGHT 600
#define WINDOW_WIDTH 1200

void InitCUDA()
{
    // We have to call cudaGLSetGLDevice if we want to use OpenGL interoperability.
	cudaGLSetGLDevice(0);
}

void Visualization(const char * input_image_name)
{
	Scene scene;
	scene.graphic_objects.reserve(4);

	Camera * camera = new Camera;
	camera->heightAngle = PI / 6.f;
	camera->aspectRatio = static_cast<float>(WINDOW_HEIGHT) / static_cast<float>(WINDOW_WIDTH);
	camera->direction = Point<3>(0.f, 0.f, -1.f);
	camera->position = Point<3>(0.f, 0.f, 6.f);
	camera->right = Point<3>(1.f, 0.f, 0.f);// right=direction x up
	camera->up = Point<3>(0.f,1.f, 0.f);
	scene.camera = camera;

	scene.graphic_objects.push_back(camera);
	Image * img = ReadImagePNG(input_image_name);	
	
	Canvas_Rectangle * reference_canvas = new Canvas_Rectangle(Point<3>(-1.f, -static_cast<double>(img->height) / static_cast<double>(img->width), 0.f), Point<3>(1.f, 0.f, 0.f),
						Point<3>(0.f, 1.f, 0.f), 2.f, 2.f*static_cast<double>(img->height) / static_cast<double>(img->width));
						reference_canvas->reference_vertex += Point<3>(1.2f, 0.f, 0.f);

	Canvas_Rectangle * roi_canvas = new Canvas_Rectangle(Point<3>(-1.f, -static_cast<double>(img->height) / static_cast<double>(img->width), 0.f), Point<3>(1.f, 0.f, 0.f),
						Point<3>(0.f, 1.f, 0.f), 2.f, 2.f*static_cast<double>(img->height) / static_cast<double>(img->width));
					    roi_canvas->reference_vertex += Point<3>(-1.2f, 0.f, 0.f);

	Array_Hierarchy * zoom_object = new Array_Hierarchy(12,img->height,img->width,img->colors,roi_canvas,reference_canvas);

	ZoomVisualization * zoom_vis = new ZoomVisualization(img->height,img->width,img->colors,zoom_object,reference_canvas);
	
	scene.graphic_objects.push_back(zoom_vis);

	scene.keyboard_active_object = 1;
	scene.mouse_active_object = 1;
	Window::View(&scene, WINDOW_WIDTH,WINDOW_HEIGHT);
}

int main(int argc, char * argv[])
{
	if (argc != 2)
	{
		printf(" Image Processing <image_name>\n");
		return 0;
	}
	else
	{
	InitCUDA();
	Visualization(argv[1]);
	return 0;
	}
}
