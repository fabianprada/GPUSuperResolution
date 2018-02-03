#include <string.h>
#include "time.h"
#include <stdlib.h>
#include "window.h"
#include "math.h"
#define PI 3.1415926

const char Window::KEY_ESCAPE = '\033';
Scene* Window::scene = NULL;
Mouse Window::mouse;
//int Window::curveFit = LINEAR;
int Window::isVisible;
int Window::frameCount;
double Window::startTime;
//double Window::radius;
double Window::frameRate;
double Window::frameCountStart;
Prompt Window::prompt;

/** This function prints out the state of the OpenGL error. */
int Window::PrintError(const int& showNoError){
	int x, y;
	int e = 1;
	switch (glGetError()){
	case GL_NO_ERROR:
		if (showNoError){ printf("No error\n"); }
		e = 0;
		break;
	case GL_INVALID_ENUM:
		printf("Invalid Enum\n");
		break;
	case GL_INVALID_VALUE:
		printf("Invalid Value\n");
		break;
	case GL_INVALID_OPERATION:
		printf("Invalid Operation\n");
		break;
	case GL_STACK_OVERFLOW:
		printf("Stack Overflow\n");
		break;
	case GL_STACK_UNDERFLOW:
		printf("Stack Underflow\n");
		break;
	case GL_OUT_OF_MEMORY:
		printf("Out of memory\n");
		break;
	}
	if (!showNoError && !e){ return 0; }

	glGetIntegerv(GL_MATRIX_MODE, &x);
	if (x == GL_MODELVIEW){
		glGetIntegerv(GL_MODELVIEW_STACK_DEPTH, &x);
		glGetIntegerv(GL_MAX_MODELVIEW_STACK_DEPTH, &y);
		printf("Modelview Matrix: %d/%d\n", x, y);
	}
	else if (x == GL_PROJECTION){
		glGetIntegerv(GL_PROJECTION_STACK_DEPTH, &x);
		glGetIntegerv(GL_MAX_PROJECTION_STACK_DEPTH, &y);
		printf("Projection Matrix: %d/%d\n", y);
	}
	return e;
}

/** This function writes out the specified string, left-aligned, at the specified location, onto the OpenGL window. */
void Window::WriteLeftString(const int& x, const int& y, const char* str){
	GLint vp[4];
	GLint dt = glIsEnabled(GL_DEPTH_TEST);
	GLint lm = glIsEnabled(GL_LIGHTING);
	GLint mm;

	glGetIntegerv(GL_VIEWPORT, vp);
	glGetIntegerv(GL_MATRIX_MODE, &mm);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, vp[2], 0, vp[3], 0, 1);

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glColor4f(0, 0, 0, 1);
	glRasterPos2f(x, y);
	int len = strlen(str);
	for (int i = 0; i<len; i++){ glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]); }
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	if (dt){ glEnable(GL_DEPTH_TEST); }
	if (lm){ glEnable(GL_LIGHTING); }
	glMatrixMode(mm);
}
/** This function writes out the specified string, right-aligned, at the specified location, onto the OpenGL window. */
void Window::WriteRightString(const int& x, const int& y, const char* str){
	GLint vp[4];
	glGetIntegerv(GL_VIEWPORT, vp);

	WriteLeftString(vp[2] - x - glutBitmapLength(GLUT_BITMAP_HELVETICA_18, (unsigned char*)str), y, str);
}

/** This function is called when no events need to be processed. */
void Window::IdleFunction(void){
	// Update the parameter values
	scene->setCurrentTime(GetTime() - startTime);
	// Just draw the scene again
	if (isVisible){ glutPostRedisplay(); }
}
/** This function is called when the visibility state of the window changes. */
void Window::VisibilityFunction(int state){
	if (state == GLUT_VISIBLE){ isVisible = 1; }
	else if (state == GLUT_NOT_VISIBLE){ isVisible = 0; }
}

//////////////////////////
// Mouse event handlers //
//////////////////////////
/** This function is called when the state of a mouse button is changed */
void Window::MouseFunction(int button, int state, int x, int y){
	mouse.update(button, state, x, y);
}
/** This function is called when one of the mouse buttons is depressed and the mouse is moved. */
void Window::MotionFunction(int x, int y){
	Point<2> d = mouse.move(x, y);

	MouseEvent mouse_event;
	mouse_event.displacement = d;
	mouse_event.position = Point<2>(mouse.endX, mouse.endY);
	if (mouse.middleDown){
		mouse_event.active_button = MIDDLE_BUTTON;
	}
	else if(mouse.rightDown){
		mouse_event.active_button = RIGHT_BUTTON;
	}
	else if (mouse.leftDown){
		mouse_event.active_button = LEFT_BUTTON;
	}
	scene->HandleMouseEvent(mouse_event);
	glutPostRedisplay();
}
/** This function is called when the mouse is moved moved but no buttons are depressed. */
void Window::PassiveMotionFunction(int x, int y){
	mouse.move(x, y);
	glutPostRedisplay();
}

/////////////////////////////
// Keyboard event handlers //
/////////////////////////////
/** This function is called when the user presses a key. */
void Window::KeyboardFunction(unsigned char c, int x, int y){
	char temp[500];
	//Image32 img;

	if (prompt.isActive)
	{
		size_t len = strlen(prompt.string_value);
		if (c == KEY_BS)
		{
			if (len>prompt.initial_length) prompt.string_value[len - 1] = 0;
		}
		else if (c == KEY_ENTER)
		{
			scene->HandleKeyboardEvent(prompt);
			prompt.isActive = false;
			prompt.string_value[0]=0;
			prompt.initial_length = 0;
			prompt.active_command = 0;
		}
		else if (c == KEY_CTRL_C)
		{
			prompt.isActive = false;
			prompt.string_value[0] = 0;
			prompt.initial_length = 0;
			prompt.active_command = 0;
		}
		else if (c >= 32 && c <= 126) // ' ' to '~'
		{
			prompt.string_value[len] = c;
			prompt.string_value[len + 1] = 0;
		}
	}
	else
	{
		prompt.active_command = c;
		scene->HandleKeyboardEvent(prompt);
	}
	glutPostRedisplay();
}
/** This function is called when the user presses one of the special keys. */
void Window::SpecialFunction(int key, int x, int y){
	switch (key){
	case GLUT_KEY_F1:
		break;
	case GLUT_KEY_F2:
		break;
	case GLUT_KEY_F3:
		break;
	case GLUT_KEY_F4:
		break;
	case GLUT_KEY_F5:
		break;
	case GLUT_KEY_F6:
		break;
	case GLUT_KEY_F7:
		break;
	case GLUT_KEY_F8:
		break;
	case GLUT_KEY_F9:
		break;
	case GLUT_KEY_F10:
		break;
	case GLUT_KEY_F11:
		break;
	case GLUT_KEY_F12:
		break;
	case GLUT_KEY_UP:
		break;
	case GLUT_KEY_DOWN:
		break;
	case GLUT_KEY_LEFT:
		break;
	case GLUT_KEY_RIGHT:
		break;
	case GLUT_KEY_PAGE_UP:
		break;
	case GLUT_KEY_PAGE_DOWN:
		break;
	case GLUT_KEY_HOME:
		break;
	case GLUT_KEY_END:
		break;
	case GLUT_KEY_INSERT:
		break;
	}
	glutPostRedisplay();
}

/////////////////////////
// Menu event handlers //
/////////////////////////
static bool drawWithOpenGL = true;
/** This function is called when the user updates the draw mode in the drop-down menu. */
void Window::DrawModeMenu(int entry){
	if (entry == GL_POINT || entry == GL_LINE || entry == GL_FILL)
	{
		drawWithOpenGL = true;
		glPolygonMode(GL_FRONT_AND_BACK, entry);
		glutPostRedisplay();
	}
	else {
		drawWithOpenGL = false;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glutPostRedisplay();
	}
}

/** This function is called when the user updates the cull mode in the drop-down menu. */
void Window::CullModeMenu(int entry){
	if (entry == NO_CULL){ glDisable(GL_CULL_FACE); }
	else{
		glEnable(GL_CULL_FACE);
		if (entry == CULL_BACK_FACE){ glCullFace(GL_BACK); }
		if (entry == CULL_FRONT_FACE){ glCullFace(GL_FRONT); }
	}
	glutPostRedisplay();
}
/** This function is called when the user updates the curve fitting method in the drop-down menu. */

/** This function is called when the user selects one of the main menu options in the drop-down menu. */
void Window::MainMenu(int entry){
	if (!entry){ exit(0); }
}

/////////////////////////

/**  This function draws the OpenGL window. */
void Window::DisplayFunction(void){
	char temp[500];
	double t;

	// Clear the color and depth buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	float d = scene->radius + (scene->camera->position - scene->center).norm();

	// Draw the perspective projection (to get good front and back clipping planes, we need to know the rough
	// size of the model)
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(scene->camera->heightAngle*180.0 / PI, scene->camera->aspectRatio, .01*d, 2 * d);
//	gluPerspective(scene->camera->heightAngle*180.0 / PI, scene->camera->aspectRatio, .01*d, 2 * d);

	// Draw the Scene
	GLint drawMode[2];
	glGetIntegerv(GL_POLYGON_MODE, drawMode);
	if (drawMode[0] == GL_FILL){ scene->drawOpenGL(); }
	else{
		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glPolygonOffset(1, 2);
		glDisable(GL_BLEND);
		scene->drawOpenGL();
		glDisable(GL_POLYGON_OFFSET_FILL);
		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		glPolygonMode(GL_FRONT_AND_BACK, drawMode[0]);
		glDisable(GL_BLEND);
		glEnable(GL_LINE_SMOOTH);
		glLineWidth(2);
		scene->drawOpenGL();
	}


	// Update the frame rate value on the tenth frame
	frameCount++;
	if (frameCount == 10){
		frameCount = 0;
		t = frameCountStart;
		frameCountStart = GetTime();
		frameRate = 10 / (frameCountStart - t);
	}
	//sprintf(temp, "%.1f fs", frameRate);
	WriteLeftString(1, 2, prompt.string_value);

	// Write out the mouse position
	//sprintf(temp, "(%3d , %3d)", mouse.endX, mouse.endY);
	sprintf(temp, "%.1f fs", frameRate);
	WriteRightString(1, 2, temp);

	// Swap the back and front buffers
	glutSwapBuffers();
}

/** This function is called when the window is resized. */
void Window::ReshapeFunction(int width, int height)
{
	GLint viewPort[4];

	glViewport(0, 0, width, height);
	glGetIntegerv(GL_VIEWPORT, viewPort);
	scene->camera->aspectRatio = (float)width / (float)height;
	glutPostRedisplay();
}

/** This function instantiates the OpenGL window, reading in the Scene from the specified file
* and setting the initial OpenGL window size. The function never returns! */
void Window::View(Scene* s, int width, int height){
	int drawMenu;
	int cullMenu;
//	int curveFitMenu;
	int argc = 1;
	char* argv[] = { "foo" };
	startTime = GetTime();

	scene = s;
	prompt.isActive = false;
	frameCountStart = GetTime();

	// Initialize the OpenGL context
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("GPU SuperResolution");

	if (glewInit() != GLEW_OK) fprintf(stderr, "glewInit failed. Exiting...\n");

	glutInitWindowSize(width, height);
	glutInitWindowPosition(0, 0);
	glClearColor(1.f,1.f,1.f,1.0);
	//glClearColor(0.f, 0.f, 0.f, 0.f);
	glutReshapeWindow(width, height);

	// Initialize the event handlers
	glutDisplayFunc(DisplayFunction);
	glutReshapeFunc(ReshapeFunction);
	glutKeyboardFunc(KeyboardFunction);
	glutSpecialFunc(SpecialFunction);
	glutMouseFunc(MouseFunction);
	glutMotionFunc(MotionFunction);
	glutPassiveMotionFunc(PassiveMotionFunction);
	glutVisibilityFunc(VisibilityFunction);
	glutIdleFunc(IdleFunction);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
	glDepthMask(GL_TRUE);
	//glDisable(GL_BLEND);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glCullFace(GL_BACK);
	glEnable(GL_CULL_FACE);
	// Set up all the necessary stuff for the RayScene
	scene->setupOpenGL();
	// Fall into the main loop
	glutMainLoop();
}
