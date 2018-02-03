#ifndef WINDOW_INCLUDE
#define WINDOW_INCLUDE

#include "graphic-object.h"
#include "mouse.h"

#define KEY_UPARROW		101
#define KEY_DOWNARROW	103
#define KEY_LEFTARROW	100
#define KEY_RIGHTARROW	102
#define KEY_PGUP		104
#define KEY_PGDN		105
#define KEY_CTRL_C  3
#define KEY_BS      8
#define KEY_ENTER  13
#define KEY_ESC    27

class Window{
public:
	/** The ascii value of the escape character */
	const static char KEY_ESCAPE;
	enum{
		/** A value indicating that no culling should be used */
		NO_CULL,
		/** A value indicating that back-face culling should be used */
		CULL_BACK_FACE,
		/** A value indicating that front-face culling should be used */
		CULL_FRONT_FACE
	};

	/** A pointer to the scene that is drawn.*/
	static Scene* scene;

	/** The current state of the mouse */
	static Mouse mouse;
	/** The center of the bounding box of the drawn scene */
	//static Point<3> center;
	///** The radius of the bounding box of the drawn scene */
	//static double radius;
	/** Is the window visible */
	static int isVisible;

	/** The number of frames that have been drawn since the last
	* updated frame rate */
	static int frameCount;
	/** The time of the last updated frame rate */
	static double frameCountStart;
	/** The current frame rate */
	static double frameRate;

	static double startTime;

	//static int curveFit;
	static Prompt prompt;

	/** This function prints out the state of the OpenGL error. */
	static int PrintError(const int& showNoError = 0);
	/** This function writes out the specified string, left-aligned, at the specified location, onto the OpenGL window. */
	static void WriteLeftString(const int& x, const int& y, const char* str);
	/** This function writes out the specified string, right-aligned, at the specified location, onto the OpenGL window. */
	static void WriteRightString(const int& x, const int& y, const char* str);
	//static void WriteString(int x, int y, void* font, bool bottom, bool left, const char* format, ...);

	/** This function reads the current frame buffer and sets the pixels of the image accordingly. */
	//static int TakeSnapshot(Image32& img);
	/** This function is called when no events need to be processed. */
	static void IdleFunction(void);
	/** This function is called when the visibility state of the window changes. */
	static void VisibilityFunction(int state);

	//////////////////////////
	// Mouse event handlers //
	//////////////////////////
	/** This function is called when the state of a mouse button is changed */
	static void MouseFunction(int button, int state, int x, int y);
	/** This function is called when one of the mouse buttons is depressed and the mouse is moved. */
	static void MotionFunction(int x, int y);
	/** This function is called when the mouse is moved moved but no buttons are depressed. */
	static void PassiveMotionFunction(int x, int y);

	/////////////////////////////
	// Keyboard event handlers //
	/////////////////////////////
	/** This function is called when the user presses a key. */
	static void KeyboardFunction(unsigned char c, int x, int y);
	/** This function is called when the user presses one of the special keys. */
	static void SpecialFunction(int key, int x, int y);

	/////////////////////////
	// Menu event handlers //
	/////////////////////////
	/** This function is called when the user updates the draw mode in the drop-down menu. */
	static void DrawModeMenu(int entry);
	/** This function is called when the user updates the cull mode in the drop-down menu. */
	static void CullModeMenu(int entry);
	/** This function is called when the user updates the curve fitting method in the drop-down menu. */
	static void CurveFitMenu(int entry);
	/** This function is called when the user selects one of the main menu options in the drop-down menu. */
	static void MainMenu(int entry);
	/**  This function draws the OpenGL window. */
	static void drawScene(float jitterX, float jitterY);
	static void DisplayFunction(void);
	/** This function is called when the window is resized. */
	static void ReshapeFunction(int width, int height);


	static void raytraceScene(void);

	/** This function instantiates the OpenGL window, reading in the RayScene from the specified file
	* and setting the initial OpenGL window size. The function never returns! */
	static void View(Scene* scene, int width, int height);
};

#endif //WINDOW_INCLUDE