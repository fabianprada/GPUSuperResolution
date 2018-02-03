#ifndef MOUSE_INCLUDED
#define MOUSE_INCLUDED

/** This class represents the current state of the mouse */

template< int dimension>
class Point;

class Mouse{
public:
	/** Was the shift button depressed when the mouse state was updated */
	int shiftDown;

	/** Was the control button depressed when the mouse state was updated */
	int ctrlDown;

	/** Was the alt button depressed when the mouse state was updated */
	int altDown;

	/** Was the left mouse button depressed when the mouse state was updated */
	int leftDown;

	/** Was the middle mouse button depressed when the mouse state was updated */
	int middleDown;

	/** Was the right mouse button depressed when the mouse state was updated */
	int rightDown;

	/** Was the scroll wheel depressed when the mouse state was updated */
	int scrollDown;

	/** The screen coordinates when the mouse state was updated */
	int startX,startY;

	/** The screen coordinates when the mouse was last moved */
	int endX,endY;

	
	/** This constructor instantiates the mouse members */
	Mouse(void);

	/** This method updates the state of the mouse when the state of a mouse button is changed */
	void update(int button,int state,int x, int y);

	/** This method udpates the state of the mouse when the mouse is moved. The returned value
	  * indicates the distance (in screen coordinates) that the mouse has moved */
	Point<2> move(int x, int y);
};
#endif // MOUSE_INCLUDED