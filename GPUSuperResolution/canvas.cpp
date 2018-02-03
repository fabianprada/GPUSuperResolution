#include "graphic-object.h"

void Canvas_Rectangle::DrawTexture(GLuint texture_buffer)
{
	glEnable(GL_TEXTURE_2D);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glBindTexture(GL_TEXTURE_2D, texture_buffer);

	glBegin(GL_QUADS);
	Point<3> vertex_position = reference_vertex;
	glTexCoord2d(0.f, 0.f);
	glVertex3d(vertex_position[0], vertex_position[1], vertex_position[2]);

	vertex_position += axis[0] * edges_lenght[0];
	glTexCoord2d(1.f, 0.f);
	glVertex3d(vertex_position[0], vertex_position[1], vertex_position[2]);

	vertex_position += axis[1] * edges_lenght[1];
	glTexCoord2d(1.f, 1.f);
	glVertex3d(vertex_position[0], vertex_position[1], vertex_position[2]);

	vertex_position -= axis[0] * edges_lenght[0];
	glTexCoord2d(0.f, 1.f);
	glVertex3d(vertex_position[0], vertex_position[1], vertex_position[2]);

	glEnd();

	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisable(GL_TEXTURE_2D);
}

void Canvas_Rectangle::DrawTriangles(GLuint & position_buffer,const unsigned int num_triangles, Point<3> p_color){

	 glMatrixMode(GL_MODELVIEW);
     glPushMatrix();
	 
	 glTranslatef(reference_vertex[0],reference_vertex[1],reference_vertex[2] + (edges_lenght[0] + edges_lenght[1]) / 10000.f);
	 glScalef(edges_lenght[0],edges_lenght[1],1.f);
	 
	glColor3f(p_color[0], p_color[1], p_color[2]);

	glEnableClientState(GL_VERTEX_ARRAY);
	glBindBuffer(GL_ARRAY_BUFFER, position_buffer);
	glVertexPointer(3, GL_FLOAT, 0, 0);
	glDrawArrays(GL_TRIANGLES, 0, 3 * num_triangles);
	glDisableClientState(GL_VERTEX_ARRAY);

	glPopMatrix();
}

void Canvas_Rectangle::DrawBox(Point<2> p_canvas_coordinates, double p_canvas_scale, Point<3> p_color, double p_width)
{
	glColor3f(p_color[0], p_color[1], p_color[2]);
	glLineWidth(p_width);
	Point<3> canvas_pos;
	Point<3> world_pos;
	glBegin(GL_LINE_LOOP);

	canvas_pos = Point<3>(p_canvas_coordinates[0], p_canvas_coordinates[1],(edges_lenght[0] + edges_lenght[1]) / 1000.f);
	world_pos = CanvasToWorld(canvas_pos[0], canvas_pos[1], canvas_pos[2]);
	glVertex3f(world_pos[0], world_pos[1], world_pos[2]);

	canvas_pos += Point<3>(p_canvas_scale, 0.f, 0.f);
	world_pos =CanvasToWorld(canvas_pos[0], canvas_pos[1], canvas_pos[2]);
	glVertex3f(world_pos[0], world_pos[1], world_pos[2]);

	canvas_pos += Point<3>(0.f, p_canvas_scale, 0.f);
	world_pos = CanvasToWorld(canvas_pos[0], canvas_pos[1], canvas_pos[2]);
	glVertex3f(world_pos[0], world_pos[1], world_pos[2]);

	canvas_pos += Point<3>(-p_canvas_scale, 0.f, 0.f);
	world_pos = CanvasToWorld(canvas_pos[0], canvas_pos[1], canvas_pos[2]);
	glVertex3f(world_pos[0], world_pos[1], world_pos[2]);

	glEnd();
}
