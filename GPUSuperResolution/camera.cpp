#include "graphic-object.h"

void Camera::drawOpenGL(Scene * scene)
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(position[0], position[1], position[2], position[0] + direction[0], position[1] + direction[1], position[2] + direction[2], up[0], up[1], up[2]);
}

class Scene;

void Camera::HandleMouseEvent(MouseEvent & mouse_event, Scene * scene)
{
	if (mouse_event.active_button == MIDDLE_BUTTON){
		moveForward(scene->radius / 2500.f * mouse_event.displacement[1]);
	}
	else if (mouse_event.active_button == RIGHT_BUTTON)
	{
		rotateUp(scene->center, 0.001f*mouse_event.displacement[0]);
		rotateRight(scene->center, 0.001f*mouse_event.displacement[1]);
	}
	else if (mouse_event.active_button == LEFT_BUTTON){
		moveRight(-scene->radius / 2500.f * mouse_event.displacement[0]);
		moveUp(scene->radius / 2500.f *  mouse_event.displacement[1]);
	}
}

void Camera::rotateUp(Point<3> center, float angle)
{
	right = CentralizedRotation(Point<3>(), up, right, angle);
	direction = CentralizedRotation(Point<3>(), up, direction, angle);
	position = CentralizedRotation(center, up, position, angle);
}

void Camera::rotateDirection(Point<3> center, float angle)
{
	right = CentralizedRotation(Point<3>(), direction, right, angle);
	up = CentralizedRotation(Point<3>(), direction, up, angle);
	position = CentralizedRotation(center, direction, position, angle);
}

void Camera::rotateRight(Point<3> center, float angle)
{
	up = CentralizedRotation(Point<3>(), right, up, angle);
	direction = CentralizedRotation(Point<3>(), right, direction, angle);
	position = CentralizedRotation(center, right, position, angle);
}
void Camera::moveForward(float dist){
	position += (direction*dist);
}
void Camera::moveRight(float dist){
	position += (right*dist);
}
void Camera::moveUp(float dist){
	position += (up*dist);
}