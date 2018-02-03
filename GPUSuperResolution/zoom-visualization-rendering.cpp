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


#include "zoom-visualization.h"

void ZoomVisualization::setupOpenGL()
{
	GLTools::CreateTexture(color_texture_buffer,width,height,color);
	zoom_object->setupOpenGL();
}


void ZoomVisualization::drawOpenGL(Scene * scene)
{
	reference_canvas->DrawTexture(color_texture_buffer);
	reference_canvas->DrawBox(roi_absolute_coordinate, roi_absolute_scale,Point<3>(1.f,1.f,0.f));
	zoom_object->drawROI(roi_absolute_coordinate,roi_absolute_scale);
}

void ZoomVisualization::HandleMouseEvent(MouseEvent  & mouse_event, Scene * scene)
{
	if (mouse_event.active_button == MIDDLE_BUTTON){
		double current_scale = roi_absolute_scale;
		double new_scale = current_scale *exp(mouse_event.displacement[1] / 500.f);
		if(new_scale > 1.f){
		   new_scale = 1.f;
		}
		roi_absolute_scale = new_scale;
		roi_absolute_coordinate -= Point<2>(1.f, 1.f)*(new_scale - current_scale) / 2.f;
	}
	else if (mouse_event.active_button == LEFT_BUTTON){
		roi_absolute_coordinate.coord[0] += mouse_event.displacement[0] / 2500.f;
		roi_absolute_coordinate.coord[1] += mouse_event.displacement[1] / 2500.f;
	}
	else if (mouse_event.active_button == RIGHT_BUTTON){
		/*camera->moveRight(mouse_event.displacement[0] / 2500.f);
		camera->moveUp(mouse_event.displacement[1] / 2500.f);*/
	}
}

void ZoomVisualization::HandleKeyboardEvent(Prompt  & prompt, Scene * scene)
{
		zoom_object->HandleKeyboardEvent(prompt,scene);
}


