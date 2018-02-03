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


#include "array-hierarchy.h"

template<int dimension>
void ZoomVisualization<dimension>::ExtractROI()
{
	roi->getValuesFromRegionUpsampling(array_hierarchy->augmented_roi->color, roi_padding_factor, roi_padding_factor, 1.f - 2 * roi_padding_factor, 1.f - 2 * roi_padding_factor, 1.f);
	roi_texture_buffer_updated = false;
}


template<int dimension>
void ZoomVisualization<dimension>::setupOpenGL()
{
	GLTools::SetupTexture(reference_texture_buffer);
	array_hierarchy->levels[0]->color->setTextureValues(reference_texture_buffer);
	GLTools::SetupTexture(roi_texture_buffer);
	roi_texture_buffer_updated = false;
}

template<int dimension>
void ZoomVisualization<dimension>::drawOpenGL(Scene * scene)
{
	// Draw Textures
	reference_canvas->DrawTexture(reference_texture_buffer);

	if (!roi_texture_buffer_updated){
		//roi->setTextureValues(roi_texture_buffer);
		int roi_corner_w = static_cast<int>(floor(static_cast<double>(array_hierarchy->augmented_roi->width)*roi_padding_factor));
		int roi_corner_h = static_cast<int>(floor(static_cast<double>(array_hierarchy->augmented_roi->height)*roi_padding_factor));
		array_hierarchy->augmented_roi->color->setTextureValues(roi_texture_buffer, roi_corner_w, roi_corner_h, width, height);
		roi_texture_buffer_updated = true;
	}
	roi_canvas->DrawTexture(roi_texture_buffer);
	
	//Draw Boxes
	drawBox(reference_canvas, roi_absolute_coordinate, roi_absolute_scale);
}

template<int dimension>
void ZoomVisualization<dimension>::drawBox(Canvas_Rectangle * p_canvas, Point<2> p_canvas_coordinates, double p_canvas_scale, Point<3> p_color = Point<3>(), double p_width = 0.01f )
{
	glColor3f(p_color[0], p_color[1], p_color[2]);
	glLineWidth(p_width);
	Point<3> canvas_pos;
	Point<3> world_pos;
	glBegin(GL_LINE_LOOP);

	canvas_pos = Point<3>(p_canvas_coordinates[0], p_canvas_coordinates[1],(p_canvas->edges_lenght[0] + p_canvas->edges_lenght[1]) / 1000.f);
	world_pos = p_canvas->CanvasToWorld(canvas_pos[0], canvas_pos[1], canvas_pos[2]);
	glVertex3f(world_pos[0], world_pos[1], world_pos[2]);

	canvas_pos += Point<3>(p_canvas_scale, 0.f, 0.f);
	world_pos = p_canvas->CanvasToWorld(canvas_pos[0], canvas_pos[1], canvas_pos[2]);
	glVertex3f(world_pos[0], world_pos[1], world_pos[2]);

	canvas_pos += Point<3>(0.f, p_canvas_scale, 0.f);
	world_pos = p_canvas->CanvasToWorld(canvas_pos[0], canvas_pos[1], canvas_pos[2]);
	glVertex3f(world_pos[0], world_pos[1], world_pos[2]);

	canvas_pos += Point<3>(-p_canvas_scale, 0.f, 0.f);
	world_pos = p_canvas->CanvasToWorld(canvas_pos[0], canvas_pos[1], canvas_pos[2]);
	glVertex3f(world_pos[0], world_pos[1], world_pos[2]);

	glEnd();

}

template<int dimension>
void ZoomVisualization<dimension>::HandleMouseEvent(MouseEvent  & mouse_event, Scene * scene = 0)
{
	if (mouse_event.active_button == MIDDLE_BUTTON){
		//printf("Exp %f \n", exp(mouse_event.displacement[1] / 50.f));
		double current_scale = roi_absolute_scale;
		double new_scale = current_scale *exp(mouse_event.displacement[1] / 500.f);
		//Point<2> new_absolute_coordinate = roi_absolute_coordinate;
		//if (new_scale < 1.f){
			roi_absolute_scale = new_scale;
			Point<2> new_absolute_coordinate = roi_absolute_coordinate - Point<2>(1.f, 1.f)*(new_scale - current_scale) / 2.f;
		//}
		//else{
		//	roi_absolute_scale = 1.f;
		//	new_absolute_coordinate = Point<2>();
		//}

		array_hierarchy->Update(new_absolute_coordinate, new_scale);
		//ExtractROI();
		roi_texture_buffer_updated = false;
	}
	else if (mouse_event.active_button == RIGHT_BUTTON){
		roi_absolute_coordinate.coord[0] += mouse_event.displacement[0] / 2500.f;
		roi_absolute_coordinate.coord[1] += mouse_event.displacement[1] / 2500.f;
		array_hierarchy->Update(roi_absolute_coordinate, roi_absolute_scale);
		roi_texture_buffer_updated = false;
		//ExtractROI();
	}
	else if (mouse_event.active_button == LEFT_BUTTON){
		/*camera->moveRight(mouse_event.displacement[0] / 2500.f);
		camera->moveUp(mouse_event.displacement[1] / 2500.f);*/
	}
}