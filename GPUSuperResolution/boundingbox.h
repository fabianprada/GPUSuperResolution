#include "graphic-object.h"

class BoundingBox : public GraphicObject{
public:
	BoundingBox()
	{
		sprintf(object_name,"BOUNDING BOX");
		is_active = true;
		reference_vertex = Point<3>(-1.f, -1.f, 1.f);
		directions[0] = Point<3>(1.f, 0.f, 0.f);
		directions[1] = Point<3>(0.f, 1.f, 0.f);
		directions[2] = Point<3>(0.f, 0.f, -1.f);
		lenghts[0] = 2.f;
		lenghts[1] = 2.f;
		lenghts[2] = 2.f;

	}
	Point<3> reference_vertex;
	Point<3> directions[3];
	double lenghts[3];
	void setupOpenGL(){}
	void drawOpenGL(Scene * scene = 0);
	void HandleKeyboardEvent(Prompt  & prompt, Scene * scene = 0);
	void HandleMouseEvent(MouseEvent & mouse_event, Scene * scene = 0);

	void rotateDirection1(float angle);
	void rotateDirection2(float angle);
	void moveX(float dist);
	void moveY(float dist);
	void moveZ(float dist);

	bool ReturnDoubleProperty(const char * property_name, double & value);
	bool ReturnPointProperty(const char * property_name, Point<3> & value);

};