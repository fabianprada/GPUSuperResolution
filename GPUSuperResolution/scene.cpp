#include "graphic-object.h"


void Material::drawOpenGL()
{
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emissive);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, specularFallOff);
}

void Scene::drawOpenGL(){

	for (int i = 0; i < graphic_objects.size(); i++){
		if (graphic_objects[i]->is_active){
		graphic_objects[i]->drawOpenGL(this);
		}
	}
}

void Scene::setupOpenGL(){
	for (int i = 0; i < graphic_objects.size(); i++){
		graphic_objects[i]->setupOpenGL();
	}
}

void Scene::setCurrentTime(double t)
{
	double velocity = 0.25f;
	t *= velocity;
	double residue = t - floor(t);
	current_time = residue;
}


void Scene::ChangeKeyboardActiveObject_CallBack(Prompt & prompt)
{
	keyboard_active_object = (keyboard_active_object + 1) % graphic_objects.size();
	sprintf(prompt.string_value, "KEYBOARD ACTIVE OBJECT : %s", graphic_objects[keyboard_active_object]->object_name);
}

void Scene::ChangeMouseActiveObject_CallBack(Prompt & prompt)
{
	mouse_active_object = (mouse_active_object + 1) % graphic_objects.size();
	sprintf(prompt.string_value, "MOUSE ACTIVE OBJECT : %s", graphic_objects[mouse_active_object]->object_name);
} 

void Scene::HandleKeyboardEvent(Prompt & prompt)
{
	bool matched_key = true;
	switch (prompt.active_command){
	case 'z':
		ChangeKeyboardActiveObject_CallBack(prompt);
		break;
	case 'x':
		ChangeMouseActiveObject_CallBack(prompt);
		break;
	default:
		matched_key = false;
	}
	if (!matched_key)
	{
		graphic_objects[keyboard_active_object]->HandleKeyboardEvent(prompt, this);
	}
}

void Scene::HandleMouseEvent(MouseEvent  & mouse_event)
{
	graphic_objects[mouse_active_object]->HandleMouseEvent(mouse_event,this);
}

bool Scene::ObjectsCommunicationDouble(const char * object_name, const char * property_name, double & value)
{
	int object_index = -1;
	int i = 0;
	while (i < graphic_objects.size() && object_index==-1){
		if (strcmp(graphic_objects[i]->object_name, object_name) == 0){
			object_index = i;
		}
		i++;
	}
	if (object_index == -1)
		return false;
	else
		return	graphic_objects[object_index]->ReturnDoubleProperty(property_name, value);
}

bool Scene::ObjectsCommunicationInt(const char * object_name, const char * property_name, int & value)
{
	int object_index = -1;
	int i = 0;
	while (i < graphic_objects.size() && object_index == -1){
		if (strcmp(graphic_objects[i]->object_name, object_name) == 0){
			object_index = i;
		}
		i++;
	}
	if (object_index == -1)
		return false;
	else
		return	graphic_objects[object_index]->ReturnIntProperty(property_name, value);
}

bool Scene::ObjectsCommunicationPoint(const char * object_name, const char * property_name, Point<3> & value)
{
	int object_index = -1;
	int i = 0;
	while (i < graphic_objects.size() && object_index == -1){
		if (strcmp(graphic_objects[i]->object_name, object_name) == 0){
			object_index = i;
		}
		i++;
	}
	if (object_index == -1)
		return false;
	else
		return	graphic_objects[object_index]->ReturnPointProperty(property_name, value);
}
