#include "graphic-object.h"

void Light::drawOpenGL(Scene * scene)
{
	GLfloat light_position[4] = { position[0], position[1], position[2], 0.f };
	if (type == DIRECTIONAL_LIGHT)
	{
		light_position[3] = 1.f;
	}
	GLfloat light_color[4] = { color[0], color[1], color[2], 1.f };
	glLightfv(GL_LIGHT0 + index, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0 + index, GL_DIFFUSE, light_color);
	glLightfv(GL_LIGHT0 + index, GL_SPECULAR, light_color);
	glEnable(GL_LIGHT0 + index);

	glPointSize(3.f);
	glEnable(GL_POINT_SMOOTH);
	glBegin(GL_POINTS);
	glColor3f(1.f - color[0], 1.f - color[1], 1.f - color[2]);
	glVertex3f(position[0], position[1], position[2]);
	glEnd();

}

void Light::rotateX(Point<3> center, float angle)
{
	position = CentralizedRotation(center, Point<3>(1.f, 0.f, 0.f), position, angle);
}

void Light::rotateY(Point<3> center, float angle)
{
	position = CentralizedRotation(center, Point<3>(0.f, 1.f, 0.f), position, angle);
}

void Light::moveX(float dist)
{
	position += Point<3>(1.f, 0.f, 0.f)*dist;
}
void Light::moveY(float dist)
{
	position += Point<3>(0.f, 1.f, 0.f)*dist;
}
void Light::moveZ(float dist)
{
	position += Point<3>(0.f, 0.f, -1.f)*dist;
}

void Light::HandleMouseEvent(MouseEvent & mouse_event, Scene * scene)
{
	if (mouse_event.active_button == MIDDLE_BUTTON){
		rotateX(scene->center, 0.001f*mouse_event.displacement[0]);
		rotateY(scene->center, 0.001f*mouse_event.displacement[1]);
	}
	else if (mouse_event.active_button == RIGHT_BUTTON)
	{
		moveZ(scene->radius / 2500.f * mouse_event.displacement[1]);
	}
	else if (mouse_event.active_button == LEFT_BUTTON){
		moveX(-scene->radius / 2500.f * mouse_event.displacement[0]);
		moveY(scene->radius / 2500.f *  mouse_event.displacement[1]);
	}
}

