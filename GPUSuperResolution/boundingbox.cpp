#include "boundingbox.h"
#include "scene.h"

void BoundingBox::drawOpenGL(Scene * scene)
{
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	for (int i = 0; i < 3; i++)
	{
		int index_A = i;
		int index_B = (i+1)%3;
		int index_C = (i+2)%3;
		Point<3> current_vertex = reference_vertex;
		glBegin(GL_QUADS);
		glColor4f(1.f, 1.f, 0.f, 0.5f);
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);
		current_vertex += directions[index_A] * lenghts[index_A];
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);
		current_vertex += directions[index_B] * lenghts[index_B];
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);
		current_vertex -= directions[index_A] * lenghts[index_A];
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);

		current_vertex += directions[index_C] * lenghts[index_C];
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);
		current_vertex += directions[index_A] * lenghts[index_A];
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);
		current_vertex -= directions[index_B] * lenghts[index_B];
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);
		current_vertex -= directions[index_A] * lenghts[index_A];
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);
		glEnd();
	}
	glDisable(GL_BLEND);

	/*glLineWidth(0.01f);
	glEnable(GL_LINE_SMOOTH);
	glColor3f(0.f, 0.f, 0.f);
	for (int i = 0; i < 3; i++)
	{
		int index_A = i;
		int index_B = (i + 1) % 3;
		int index_C = (i + 2) % 3;

		Point<3> current_vertex = reference_vertex;
		glBegin(GL_LINE_LOOP);
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);
		current_vertex += directions[index_A] * lenghts[index_A];
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);
		current_vertex += directions[index_B] * lenghts[index_B];
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);
		current_vertex -= directions[index_A] * lenghts[index_A];
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);
		current_vertex += directions[index_C] * lenghts[index_C];
		glEnd();

		glBegin(GL_LINE_LOOP);
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);
		current_vertex += directions[index_A] * lenghts[index_A];
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);
		current_vertex -= directions[index_B] * lenghts[index_B];
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);
		current_vertex -= directions[index_A] * lenghts[index_A];
		glVertex3d(current_vertex[0], current_vertex[1], current_vertex[2]);
		glEnd();
	}*/
}

void BoundingBox::moveX(float dist)
{
	reference_vertex += Point<3>(1.f, 0.f, 0.f)*dist;
}
void BoundingBox::moveY(float dist)
{
	reference_vertex += Point<3>(0.f, 1.f, 0.f)*dist;
}
void BoundingBox::moveZ(float dist)
{
	reference_vertex += Point<3>(0.f, 1.f, 0.f)*dist;
}


void BoundingBox::rotateDirection1(float angle)
{
	directions[1] = CentralizedRotation(Point<3>(), directions[0], directions[1], angle);
	directions[2] = CentralizedRotation(Point<3>(), directions[0], directions[2], angle);
}

void BoundingBox::rotateDirection2(float angle)
{
	directions[0] = CentralizedRotation(Point<3>(), directions[1], directions[0], angle);
	directions[2] = CentralizedRotation(Point<3>(), directions[1], directions[2], angle);
}

void BoundingBox::HandleMouseEvent(MouseEvent & mouse_event, Scene * scene)
{
	if (mouse_event.active_button == MIDDLE_BUTTON){
		rotateDirection1(0.001f*mouse_event.displacement[0]);
		rotateDirection2(0.001f*mouse_event.displacement[1]);
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

void BoundingBox::HandleKeyboardEvent(Prompt & prompt, Scene * scene)
{
	switch (prompt.active_command){
	case 'q':
		lenghts[0] = lenghts[0] + 0.001f;
		break;
	case 'w':
		lenghts[0] = lenghts[0] - 0.001f > 0.f ? lenghts[0] - 0.001f : 0.f;
		break;
	case 'e':
		lenghts[1] = lenghts[1] + 0.001f;
		break;
	case 'r':
		lenghts[1] = lenghts[1] - 0.001f > 0.f ? lenghts[1] - 0.001f : 0.f;
		break;
	case 'a':
		lenghts[2] = lenghts[2] + 0.001f;
		break;
	case 's':
		lenghts[2] = lenghts[2] - 0.001f > 0.f ? lenghts[2] - 0.001f : 0.f;
		break;
	case 'v':
		is_active = false;
		break;
	}
}

bool BoundingBox::ReturnDoubleProperty(const char * property_name, double & value)
{
	if (strcmp(property_name, "lenghts[0]") == 0){
		value = lenghts[0];
		return true;
	}
	else if (strcmp(property_name, "lenghts[1]") == 0){
		value = lenghts[1];
		return true;
	}
	else if(strcmp(property_name, "lenghts[2]") == 0){
		value = lenghts[2];
		return true;
	}
	return false;
}
bool BoundingBox::ReturnPointProperty(const char * property_name, Point<3> & value)
{
	if (strcmp(property_name, "reference_vertex") == 0){
		value = reference_vertex;
		return true;
	}
	else if (strcmp(property_name, "directions[0]") == 0){
		value = directions[0];
		return true;
	}
	else if (strcmp(property_name, "directions[1]") == 0){
		value = directions[1];
		return true;
	}
	else if (strcmp(property_name, "directions[2]") == 0){
		value = directions[2];
		return true;
	}
	else return false;
}

bool InsideBox(Point<3> position, Point<3> reference_vertex, Point<3> direction0, Point<3> direction1, Point<3> direction2, double lenght0, double lenght1, double lenght2)
{
	Point<3> centralized_position = position - reference_vertex;
	double projection = centralized_position.dot(direction0);
	if (projection <0.f || projection > lenght0){
		return false;
	}
	projection = centralized_position.dot(direction1);
	if (projection <0.f || projection > lenght1){
		return false;
	}
	projection = centralized_position.dot(direction2);
	if (projection <0.f || projection > lenght2){
		return false;
	}
	return true;
}