#include "zoom.h"

void ZoomBox::drawOpenGL(Scene * scene)
{
	glColor3f(0.f, 0.f, 0.f);
	glLineWidth(0.01);
	Point<3> canvas_pos;
	Point<3> world_pos;
	glBegin(GL_LINE_LOOP);
	
	canvas_pos = canvas_coordinate;
	world_pos = reference_canvas->CanvasToWorld(canvas_pos[0], canvas_pos[1], canvas_pos[2]);
	glVertex3f(world_pos[0], world_pos[1], world_pos[2]);

	canvas_pos += Point<3>(1.f/zoom_factor, 0.f, 0.f);
	world_pos = reference_canvas->CanvasToWorld(canvas_pos[0], canvas_pos[1], canvas_pos[2]);
	glVertex3f(world_pos[0], world_pos[1], world_pos[2]);

	canvas_pos += Point<3>(0.f, 1.f / zoom_factor, 0.f);
	world_pos = reference_canvas->CanvasToWorld(canvas_pos[0], canvas_pos[1], canvas_pos[2]);
	glVertex3f(world_pos[0], world_pos[1], world_pos[2]);

	canvas_pos += Point<3>(-1.f / zoom_factor, 0.f, 0.f);
	world_pos = reference_canvas->CanvasToWorld(canvas_pos[0], canvas_pos[1], canvas_pos[2]);
	glVertex3f(world_pos[0], world_pos[1], world_pos[2]);

	glEnd();
}

void ZoomBox::moveX(float dist)
{
	canvas_coordinate += Point<3>(1.f, 0.f, 0.f)*dist;
}
void ZoomBox::moveY(float dist)
{
	canvas_coordinate += Point<3>(0.f, -1.f, 0.f)*dist;
}

void ZoomBox::HandleMouseEvent(MouseEvent & mouse_event, Scene * scene)
{
	if (mouse_event.active_button == LEFT_BUTTON){
		moveX(mouse_event.displacement[0] /2500.f);
		moveY(mouse_event.displacement[1] /2500.f);
	}
}