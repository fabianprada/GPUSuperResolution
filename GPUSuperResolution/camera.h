#include "graphic-object.h"

class Camera : public GraphicObject {
public:
	Camera(){
		sprintf(object_name, "CAMERA");
		is_active = true;
	}
	double heightAngle;
	double aspectRatio;
	Point<3> position;
	Point<3> direction;
	Point<3> up;
	Point<3> right;

	void drawOpenGL(Scene * scene = 0);
	void setupOpenGL(){}

	void HandleKeyboardEvent(Prompt  & prompt, Scene * scene = 0){}
	void HandleMouseEvent(MouseEvent & mouse_event, Scene * scene = 0);

	/** This call rotates the camera, about the axis which is parallel to the up direction of the camera, and passes
	* through the specified point. */
	void rotateUp(Point<3> center, float angle);
	/** This call rotates the camera, about the axis which is parallel to the right direction of the camera, and passes
	* through the specified point. */
	void rotateRight(Point<3> center, float angle);
	void rotateDirection(Point<3> center, float angle);
	/** This call moves the camera in the forward direction by the specified distance.*/
	void moveForward(float dist);
	/** This call moves the camera in the right direction by the specified distance.*/
	void moveRight(float dist);
	/** This call moves the camera in the up direction by the specified distance.*/
	void moveUp(float dist);
};