#include "mouse.h"
#include "geometry.h"
#include "GL\glut.h"

Mouse::Mouse(void){
	startX=endX=-1;
	startY=endY=-1;
	scrollDown=leftDown=middleDown=rightDown=0;
	shiftDown=ctrlDown=leftDown=0;
}
void Mouse::update(int button,int state,int x, int y){
	int c;

	if(state==GLUT_DOWN){
		startX=endX=x;
		startY=endY=y;
		if(button==GLUT_LEFT_BUTTON		){leftDown=1;}
		if(button==GLUT_MIDDLE_BUTTON	){middleDown=1;}
		if(button==GLUT_RIGHT_BUTTON	){rightDown=1;}
#ifdef GLUT_SCROLL_WHEEL
		if(button==GLUT_SCROLL_WHEEL	){scrollDown=1;}
#endif

		c=glutGetModifiers();
		if(c&GLUT_ACTIVE_SHIFT){shiftDown=1;}
		else{shiftDown=0;}
		if(c&GLUT_ACTIVE_CTRL){ctrlDown=1;}
		else{ctrlDown=0;}
		if(c&GLUT_ACTIVE_ALT){altDown=1;}
		else{altDown=0;}
	}
	else if(state==GLUT_UP){
		endX=x;
		endY=y;
		if(button==GLUT_LEFT_BUTTON		){leftDown=0;}
		if(button==GLUT_MIDDLE_BUTTON	){middleDown=0;}
		if(button==GLUT_RIGHT_BUTTON	){rightDown=0;}
#ifdef GLUT_SCROLL_WHEEL
		if(button==GLUT_SCROLL_WHEEL	){scrollDown=0;}
#endif
	}
}
Point<2> Mouse::move(int x,int y){
	Point<2> d;

	d.coord[0]=x-endX;
	d.coord[1]=y-endY;
	endX=x;
	endY=y;
	return d;
}
