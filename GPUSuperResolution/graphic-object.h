#ifndef GRAPHIC_OBJECT_INCLUDE
#define GRAPHIC_OBJECT_INCLUDE

#include "geometry.h"
#include "gltools.h"
#include "casting.h"


class Prompt{
public:
	Prompt()
	{
		isActive = false;
		string_value[0] = 0;
		initial_length = 0;
		active_command = 0;
	}
	bool isActive;
	char string_value[500];
	int initial_length;
	unsigned char active_command;
};

enum ActiveButton{ RIGHT_BUTTON, LEFT_BUTTON, MIDDLE_BUTTON};

class MouseEvent{
public:
	ActiveButton active_button;
	Point<2> position;
	Point<2> displacement;
};


class Material {
public:
	GLfloat ambient[3];
	GLfloat diffuse[3];
	GLfloat specular[3];
	GLfloat emissive[3];
	float specularFallOff;
	void drawOpenGL();
};

class Scene;

class GraphicObject {
public:
	bool is_active;
	virtual void setupOpenGL() = 0; // virtual () = 0; -> makes the class abstract and force function declaration in any non-abstract descendant class. 
	virtual void drawOpenGL(Scene * scene = 0) = 0;
	virtual void HandleKeyboardEvent(Prompt  & prompt, Scene * scene = 0) {}
	virtual void HandleMouseEvent(MouseEvent  & mouse_event, Scene * scene = 0) {}

	virtual bool ReturnDoubleProperty(const char * property_name, double & value){ return false; }
	virtual bool ReturnIntProperty(const char * property_name, int & value){ return false; }
	virtual bool ReturnPointProperty(const char * property_name, Point<3> & value){ return false; }
	
	virtual bool AskDoubleProperty(Scene * scene, const char * object_name, const char * property_name, double & value){ return false; }
	virtual bool AskIntProperty(Scene * scene, const char * object_name, const char * property_name, double & value){ return false; }
	virtual bool AskPointProperty(Scene * scene, const char * object_name, const char * property_name, double & value){ return false; }

	char object_name[100];
};

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

	void rotateUp(Point<3> center, float angle);
	void rotateRight(Point<3> center, float angle);
	void rotateDirection(Point<3> center, float angle);
	void moveForward(float dist);
	void moveRight(float dist);
	void moveUp(float dist);
};

enum LightType{ POINT_LIGHT, DIRECTIONAL_LIGHT };

class  Light : public GraphicObject{
public:
	Light()
	{
		sprintf(object_name, "LIGHT");
		is_active = true;
	}
	int index;
	LightType type;
	Point<3> color;
	Point<3> position;
	void drawOpenGL(Scene * scene = 0);
	void setupOpenGL(){}
	void HandleKeyboardEvent(Prompt  & prompt, Scene * scene = 0){}
	void HandleMouseEvent(MouseEvent & mouse_event, Scene * scene = 0);
	void rotateX(Point<3> center, float angle);
	void rotateY(Point<3> center, float angle);

	void moveX(float dist);
	void moveY(float dist);
	void moveZ(float dist);
};


typedef std::vector<GraphicObject *> GraphicObjects;

class Scene{
public:
	Scene()
	{
		mouse_active_object = 0;
		keyboard_active_object = 0;
		radius = 5.f;
		center = Point<3>();
	}

	double current_time;
	double radius;
	Point<3> center;

	Camera* camera;
	//Light** lights;
	//int lightNum;
	GraphicObjects graphic_objects;
	int keyboard_active_object;
	int mouse_active_object;

	//MeshObject * mesh;
	/** The root of the scene-graph */
	//class StaticRayGroup* group;
	//~Scene(void);

	/** This method returns the material with the specified index.*/
	//Material* getMaterial(int index);

	//////////////////
	// OpenGL stuff //
	//////////////////

	//void RayView(const int& width, const int& height, const int& cplx, const int& factor = RayKeyData::MATRIX);

	/** This method calls the OpenGL commands for drawing the scene. */
	void drawOpenGL();

	/** This method class the OpenGL commands to set up everything that needs to be set up prior to rendering */
	void setupOpenGL();

	/** This method updates the current time, changing the parameter values as needed */
	void setCurrentTime(double t);

	void ChangeKeyboardActiveObject_CallBack(Prompt & prompt);
	void ChangeMouseActiveObject_CallBack(Prompt & prompt);
	void HandleKeyboardEvent(Prompt  & prompt);
	void HandleMouseEvent(MouseEvent  & mouse_event);
	bool ObjectsCommunicationDouble(const char * object_name, const char * property_name, double & value);
	bool ObjectsCommunicationInt(const char * object_name, const char * property_name, int & value);
	bool ObjectsCommunicationPoint(const char * object_name, const char * property_name, Point<3> & value);
};

class Canvas_Rectangle{
public:
	Canvas_Rectangle(Point<3> p_reference_vertex, Point<3> p_axis_0, Point<3> p_axis_1, double p_edge_lenght_0, double p_edge_lenght_1){
		reference_vertex = p_reference_vertex;
		axis[0] = p_axis_0;
		axis[1] = p_axis_1;
		axis[2] = axis[0].cross(axis[1]);
		edges_lenght[0] = p_edge_lenght_0;
		edges_lenght[1] = p_edge_lenght_1;
	}
	Point<3> reference_vertex;
	Point<3> axis[3]; // right, up, front
	Point<3> CanvasToWorld(double u, double v, double w){
		return (reference_vertex + (axis[0] * edges_lenght[0] * u) + (axis[1] * edges_lenght[1] * v) + axis[2] * w);
	}
	void DrawTriangle();
	void DrawTriangles(GLuint & position_buffer,const unsigned int num_triangles, Point<3> p_color = Point<3>());
	void DrawTexture(GLuint texture_buffer);
	void DrawBox(Point<2> p_canvas_coordinates, double p_canvas_scale, Point<3> p_color = Point<3>(), double p_width = 0.01f );
	double edges_lenght[2];
};



#endif //GRAPHIC_OBJECT_INCLUDE