#include "array-zoom-in.h"

template<int dimension>
void Array_ZoomIn<dimension>::ChangeVisualizationArray_CallBack(Prompt & prompt)
{
	if (visualization_mode == INITIAL_VALUE_VISUALIZATION)
	{
		visualization_mode = CURRENT_VALUE_VISUALIZATION;
		sprintf(prompt.string_value, "VISUALIZATION MODE : CURRENT_VALUE_VISUALIZATION");
	}
	else if (visualization_mode == CURRENT_VALUE_VISUALIZATION)
	{
		visualization_mode = INITIAL_VALUE_VISUALIZATION;
		sprintf(prompt.string_value, "VISUALIZATION MODE : INITIAL_VALUE_VISUALIZATION");
	}
	/*else if (visualization_mode == GRADIENT_NORM_VISUALIZATION)
	{
		visualization_mode = HESSIAN_NORM_VISUALIZATION;
		sprintf(prompt.string_value, "VISUALIZATION MODE : HESSIAN_NORM_VISUALIZATION");
	}
	else if (visualization_mode == HESSIAN_NORM_VISUALIZATION)
	{
		visualization_mode = LAPLACIAN_VISUALIZATION;
		sprintf(prompt.string_value, "VISUALIZATION MODE :LAPLACIAN_VISUALIZATION");
	}
	else if (visualization_mode == LAPLACIAN_VISUALIZATION)
	{
		visualization_mode = INITIAL_VALUE_VISUALIZATION;
		sprintf(prompt.string_value, "VISUALIZATION MODE : INITIAL_VALUE_VISUALIZATION");
	}*/
	else
	{
		sprintf(prompt.string_value, "Undefined visualization mode");
	}
	pixel_texture_buffer_updated = false;
}

template<int dimension>
void Array_ZoomIn<dimension>::ChangeVectorField_CallBack(Prompt & prompt)
{
	if (vector_field_mode == NONE_VECTOR_FIELD)
	{
		vector_field_mode = MEAN_GRADIENT_FIELD;
		sprintf(prompt.string_value, "VECTOR FIELD MODE : MEAN_GRADIENT_FIELD");
	}
	else if (vector_field_mode == MEAN_GRADIENT_FIELD)
	{
		vector_field_mode = ADVECTION_VECTOR_FIELD;
		sprintf(prompt.string_value, "VECTOR FIELD MODE : ADVECTION_VECTOR_FIELD");
	}
	else if (vector_field_mode == ADVECTION_VECTOR_FIELD)
	{
		vector_field_mode = ACCUMULATION_VECTOR_FIELD;
		sprintf(prompt.string_value, "VECTOR FIELD MODE : ACCUMULATION_VECTOR_FIELD");
	}
	else if (vector_field_mode == ACCUMULATION_VECTOR_FIELD)
	{
		vector_field_mode = PROGRESSIVE_ADVECTION_DISPLACEMENT;
		sprintf(prompt.string_value, "VECTOR FIELD MODE : PROGRESSIVE_ADVECTION_DISPLACEMENT");
	}
	//else if (vector_field_mode == ADVECTION_VECTOR_FIELD)
	//{
	//	vector_field_mode = HESSIAN_PRINCIPAL_FIELD;
	//	sprintf(prompt.string_value, "VECTOR FIELD MODE : HESSIAN_PRINCIPAL_FIELD");
	//}
	//else if (vector_field_mode == HESSIAN_PRINCIPAL_FIELD)
	//{
	//	vector_field_mode = HESSIAN_SECONDARY_FIELD;
	//	sprintf(prompt.string_value, "VECTOR FIELD MODE : HESSIAN_SECONDARY_FIELD");
	//}
	else if (vector_field_mode == PROGRESSIVE_ADVECTION_DISPLACEMENT)
	{
		vector_field_mode = NONE_VECTOR_FIELD;
		sprintf(prompt.string_value, "VECTOR FIELD MODE : NONE_VECTOR_FIELD");
	}
	else
	{
		sprintf(prompt.string_value, "Undefined visualization mode");
	}
	vector_field_buffer_updated = false;
}


template<int dimension>
void Array_ZoomIn<dimension>::Advection_CallBack(Prompt & prompt)
{
	Advect();
	vector_field_mode = ACCUMULATION_VECTOR_FIELD;
	vector_field_buffer_updated = false;
	sprintf(prompt.string_value, "VECTOR FIELD MODE : ACCUMULATION_VECTOR_FIELD");
}

template<int dimension>
void Array_ZoomIn<dimension>::IntegrateGradientField_CallBack(Prompt & prompt)
{
	//SolveLinearSystem(GRADIENT_FIELD_LEAST_SQUARES);
	SolveFromCholesky(GRADIENT_FIELD_LEAST_SQUARES);
	vector_field_mode = NONE_VECTOR_FIELD;
	vector_field_buffer_updated = false;
	visualization_mode = CURRENT_VALUE_VISUALIZATION;
	pixel_texture_buffer_updated = false;
	sprintf(prompt.string_value, "INTEGRATED GRADIENT FIELD");
}

template<int dimension>
void Array_ZoomIn<dimension>::Advance_CallBack(Prompt & prompt)
{
	Advance();
	vector_field_mode = NONE_VECTOR_FIELD;
	vector_field_buffer_updated = false;
	visualization_mode = CURRENT_VALUE_VISUALIZATION;
	pixel_texture_buffer_updated = false;
	sprintf(prompt.string_value, "INTEGRATED GRADIENT FIELD");
}

template<int dimension>
void Array_ZoomIn<dimension>::HandleKeyboardEvent(Prompt  & prompt, Scene * scene = 0)
{
	switch (prompt.active_command){
	case 'e':
		ExtractCurrentState_CallBack(prompt);
		break;
	case 'c':
		ChangeVisualizationArray_CallBack(prompt);
		break;
	case 'v':
		ChangeVectorField_CallBack(prompt);
		break;
	case 'a':
		Advection_CallBack(prompt);
		break;
	case 'd':
		Advance_CallBack(prompt);
		break;
	case'i':
		IntegrateGradientField_CallBack(prompt);
		break;
	}
}



template<int dimension>
void Array_ZoomIn<dimension>::ExtractCurrentState_CallBack(Prompt  & prompt){
	
	ExtractCurrentState();
	pixel_texture_buffer_updated = false;
	vector_field_buffer_updated = false;
	sprintf(prompt.string_value, "EXTRACTED CURRENT FIELD");
}

template<int dimension>
void Array_ZoomIn<dimension>::setupOpenGL()
{
	setupInitialTexture(initial_texture_buffer);
	setupTexture(pixel_texture_buffer);
	setupVectorField(vector_field_buffer);
}

template<int dimension>
void Array_ZoomIn<dimension>::setupVectorField(GLuint & vector_buffer)
{
	if (!glIsBuffer(vector_buffer)){
		glGenBuffers(1, &vector_buffer);
		glBindBuffer(GL_ARRAY_BUFFER, vector_buffer);
		glBufferData(GL_ARRAY_BUFFER, 9 * array_size*sizeof(double), NULL, GL_DYNAMIC_DRAW);
	}
	vector_field_buffer_updated = false;
}


template<int dimension>
void Array_ZoomIn<dimension>::setupTexture(GLuint & texture_buffer)
{
	if (!glIsTexture(texture_buffer))
	{
		glGenTextures(1, &texture_buffer);
		glBindTexture(GL_TEXTURE_2D, texture_buffer);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	}
	pixel_texture_buffer_updated = false;
}

template<int dimension>
void Array_ZoomIn<dimension>::setupInitialTexture(GLuint & texture_buffer)
{
	if (!glIsTexture(texture_buffer))
	{
		glGenTextures(1, &texture_buffer);
		glBindTexture(GL_TEXTURE_2D, texture_buffer);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	}

	initial_field->setTextureValues(texture_buffer);
}

template<int dimension>
void Array_ZoomIn<dimension>::setTextureValues(GLuint & texture_buffer)
{
	if (visualization_mode == INITIAL_VALUE_VISUALIZATION)
	{
		initial_field->setTextureValues(texture_buffer);
	}
	else if (visualization_mode == CURRENT_VALUE_VISUALIZATION)
	{
		current_field->setTextureValues(texture_buffer);
	}
	/*else if (visualization_mode == LAPLACIAN_VISUALIZATION)
	{
		laplacian_field->setTextureValues(texture_buffer);
	}
	else if (visualization_mode == GRADIENT_NORM_VISUALIZATION)
	{
		Array2D<dimension> * gradient_norm = generateGradientNormField();
		gradient_norm->setTextureValues(texture_buffer);
		delete gradient_norm;
	}
	else if (visualization_mode == HESSIAN_NORM_VISUALIZATION)
	{
		Array2D<dimension> * hessian_norm = generateHessianNormField();
		hessian_norm->setTextureValues(texture_buffer);
		delete hessian_norm;
	}*/

	pixel_texture_buffer_updated = true;
}

template<int dimension>
void Array_ZoomIn<dimension>::setVectorFieldValues(GLuint & vector_buffer)
{
	double * coordinates = new  double[9 * array_size];

	if (vector_field_mode != PROGRESSIVE_ADVECTION_DISPLACEMENT)
	{
		Array2D<2> * vector_field = generateVectorField(vector_field_mode);
		vector_field->statistics->verbose = true;
		vector_field->statistics->setArrayStatisticsL2();

		Point<2> normalization_vector(1.f / static_cast<double>(width), 1.f / static_cast<double>(height));
		//Point<2> normalization_vector(0.5f,0.5f);
		double norm_clamping_value = vector_field->statistics->value_l2_norm_percentiles[vector_field->statistics->active_percentile];
		printf("Clamping Value %g \n", norm_clamping_value);
#pragma omp parallel for
		for (int j = 0; j < height; j++){
			for (int i = 0; i < width; i++){
				Point<2> sampling_pos = vector_field->arrayIndexTextureCoordinate(i, j);
				Point<2> vector_value = vector_field->value(i, j);
				double scalar_factor = vector_value.norm()  >  norm_clamping_value ? vector_value.norm() : norm_clamping_value;
				vector_value *= (vector_visualization_scale_factor / scalar_factor);
				Point<2> rotated_vector = Point<2>(-vector_value[1], vector_value[0]);
				Point<2> p[3];
				p[0] = sampling_pos + vector_value;
				p[1] = sampling_pos + rotated_vector*0.25f;
				p[2] = sampling_pos - rotated_vector*0.25f;
				int pixel_index = i + j*width;
				for (int k = 0; k < 3; k++){
					p[k] *= normalization_vector;
					Point<3> canvas_pos = canvas->CanvasToWorld(p[k][0], p[k][1], 0.001f);
					coordinates[9 * pixel_index + 3 * k] = canvas_pos[0];
					coordinates[9 * pixel_index + 3 * k + 1] = canvas_pos[1];
					coordinates[9 * pixel_index + 3 * k + 2] = canvas_pos[2];
				}
			}
		}
		delete vector_field;
	}
	else
	{
		Array2D<2> * displacement_field = ProgressiveAdvection();
		Point<2> normalization_vector(1.f / static_cast<double>(width), 1.f / static_cast<double>(height));
#pragma omp parallel for
		for (int j = 0; j < height; j++){
			for (int i = 0; i < width; i++){
				Point<2> sampling_pos = displacement_field->arrayIndexTextureCoordinate(i, j);
				Point<2> final_pos = displacement_field->value(i, j);
				Point<2> vector_value = final_pos - sampling_pos;
				Point<2> rotated_vector = Point<2>(-vector_value[1], vector_value[0]);
				Point<2> p[3];
				p[0] = final_pos;
				p[1] = sampling_pos + rotated_vector*0.25f;
				p[2] = sampling_pos - rotated_vector*0.25f;
				int pixel_index = i + j*width;
				for (int k = 0; k < 3; k++){
					p[k] *= normalization_vector;
					Point<3> canvas_pos = canvas->CanvasToWorld(p[k][0], p[k][1], 0.001f);
					coordinates[9 * pixel_index + 3 * k] = canvas_pos[0];
					coordinates[9 * pixel_index + 3 * k + 1] = canvas_pos[1];
					coordinates[9 * pixel_index + 3 * k + 2] = canvas_pos[2];
				}
			}
		}
		delete displacement_field;
	}

	glBindBuffer(GL_ARRAY_BUFFER, vector_buffer);
	glBufferData(GL_ARRAY_BUFFER, 9 * array_size*sizeof(double), coordinates, GL_DYNAMIC_DRAW);
	delete coordinates;
	

	vector_field_buffer_updated = true;
}


template<int dimension>
void Array_ZoomIn<dimension>::HandleMouseEvent(MouseEvent  & mouse_event, Scene * scene = 0)
{
	if (mouse_event.active_button == MIDDLE_BUTTON){
		//printf("Exp %f \n", exp(mouse_event.displacement[1] / 50.f));
		zoom_box->zoom_factor *= exp(-mouse_event.displacement[1] / 500.f);
		if (zoom_box->zoom_factor < 1.f){
			zoom_box->zoom_factor = 1.f;
		}
		SampleZoomRegion();
	}
	else if (mouse_event.active_button == RIGHT_BUTTON){
		zoom_box->moveX(mouse_event.displacement[0] / 2500.f);
		zoom_box->moveY(mouse_event.displacement[1] / 2500.f);
		SampleZoomRegion();

	}
	else if (mouse_event.active_button == LEFT_BUTTON){
		camera->moveRight(mouse_event.displacement[0] / 2500.f);
		camera->moveUp(mouse_event.displacement[1] / 2500.f);
	}
}

template<int dimension>
void Array_ZoomIn<dimension>::drawOpenGL(Scene * scene)
{
	// Draw Texture
	if (!pixel_texture_buffer_updated){
		setTextureValues(pixel_texture_buffer);
	}

	canvas->DrawTexture(pixel_texture_buffer);
	initial_texture_canvas->DrawTexture(initial_texture_buffer);
	

	//// Draw Vector Field
	if (vector_field_mode != NONE_VECTOR_FIELD){
		if (!vector_field_buffer_updated){
			setVectorFieldValues(vector_field_buffer);
		}
		glEnableClientState(GL_VERTEX_ARRAY);
		glBindBuffer(GL_ARRAY_BUFFER, vector_field_buffer);
		glVertexPointer(3, GL_DOUBLE, 0, 0);
		glColor3f(0.f, 0.f, 0.f);
		glDrawArrays(GL_TRIANGLES, 0, 3 * array_size);
		glDisableClientState(GL_VERTEX_ARRAY);
	}

	////Draw Paths 
	//if (path_is_visible){
	//	drawAllPaths();
	//}

	zoom_box->drawOpenGL();
}