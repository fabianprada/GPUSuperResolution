#ifndef ZOOMBOX_INCLUDE
#define ZOOMBOX_INCLUDE
#include "graphic-object.h"


class ZoomBox : public GraphicObject{
public:
	ZoomBox(Canvas_Rectangle * p_reference_canvas, Point<3> p_canvas_coordinate = Point<3>(0.f, 0.f, 0.f), double p_zoom_factor = 2.f){
		sprintf(object_name, "ZOOM BOX");
		is_active = true;
		
		reference_canvas = p_reference_canvas;
		canvas_coordinate = p_canvas_coordinate;
		zoom_factor = p_zoom_factor;
	}

	void setupOpenGL(){};
	void drawOpenGL(Scene * scene = 0);
	//void HandleKeyboardEvent(Prompt  & prompt, Scene * scene = 0); 
	void HandleMouseEvent(MouseEvent  & mouse_event, Scene * scene = 0);

	void moveX(float dist);
	void moveY(float dist);

	Canvas_Rectangle * reference_canvas;
	Point<3> canvas_coordinate;
	double zoom_factor;
};

#endif //ZOOMBOX_INCLUDE