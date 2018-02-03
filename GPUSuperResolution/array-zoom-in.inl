#include "array-zoom-in.h"

template<int dimension>
void Array_ZoomIn<dimension>::SampleGradientAccumulationFromCurrentField()
{
	//gradient_accumulation_w->getValuesFromDerivativeFiltering(current_field, DX);
	//gradient_accumulation_h->getValuesFromDerivativeFiltering(current_field, DY);

	printf("SAMPLING GRADIENT FIELD FROM CURRENT \n");
	gradient_accumulation_w->getValuesFromDiscreteFiltering(current_field, derivative_w);
	//gradient_accumulation_w->statistics->updateStatistics();
	gradient_accumulation_h->getValuesFromDiscreteFiltering(current_field, derivative_h);
	//gradient_accumulation_h->statistics->updateStatistics();
	printf("DONE!! \n");
}

template<int dimension>
Point<2> Array_ZoomIn<dimension>::MeanGradient(Point<2> sample_position)
{
	return Point<2>((mean_initial_field->applyDerivativeFilter(sample_position, DX))[0], (mean_initial_field->applyDerivativeFilter(sample_position, DY))[0]);
}

template<int dimension>
Matrix2D Array_ZoomIn<dimension>::MeanHessian(Point<2> sample_position)
{
	Matrix2D matrix;
	matrix.m[0][0] = (mean_initial_field->applyDerivativeFilter(sample_position, DXDX))[0];
	matrix.m[0][1] =  matrix.m[1][0] = (mean_initial_field->applyDerivativeFilter(sample_position, DXDY))[0];
	matrix.m[1][1] = (mean_initial_field->applyDerivativeFilter(sample_position, DYDY))[0];
	return matrix;
}


template<int dimension>
Point<2> Array_ZoomIn<dimension>::AdvectionVector(Point<2> sample_position)
{
	if (advection_vector_mode == MEAN_GRADIENT_PROJECTED_GRADIENT || advection_vector_mode == GRADIENT_DIRECTION_SECOND_DERIVATIVE || advection_vector_mode == NORMALIZED_MEAN_GRADIENT_PROJECTED_GRADIENT){
		//Mean Gradient
		Point<2> mean_gradient = MeanGradient(sample_position);
		Matrix2D mean_hessian_matrix = MeanHessian(sample_position);

		double mean_gradient_norm = mean_gradient.norm();

		if (mean_gradient_norm > 0.f)
		{
			double direction = (mean_hessian_matrix*mean_gradient).dot(mean_gradient);
			if (advection_vector_mode == MEAN_GRADIENT_PROJECTED_GRADIENT){
				return mean_gradient*(direction / (mean_gradient_norm*mean_gradient_norm));
			}
			else if (advection_vector_mode == GRADIENT_DIRECTION_SECOND_DERIVATIVE){
				return mean_gradient*(direction / (mean_gradient_norm*mean_gradient_norm*mean_gradient_norm));
			}
			else if (advection_vector_mode == NORMALIZED_MEAN_GRADIENT_PROJECTED_GRADIENT){
				direction = direction > 0.f ? 1.f : -1.f;
				return mean_gradient*(direction / mean_gradient_norm);
			}
		}
		else{
			return Point<2>();
		}
	}
	else{
		printf("Unknown Advection Mode \n");
		return Point<2>();
	}
}

//template<int dimension>
//Array2D<2> * Array_ZoomIn<dimension>::generateDisplacementField(VectorFieldMode vector_field)
//{
//	if (vector_field == PROGRESSIVE_ADVECTION_DISPLACEMENT){
//		Array2D<2> * vector_field_array = ProgessiveAdvection();
//#pragma omp parallel for
//		for (int j = 0; j < height; j++){
//			for (int i = 0; i < width; i++){
//				vector_field_array->value(i, j) -= vector_field_array->arrayIndexTextureCoordinate(i, j);
//			}
//		}
//		return vector_field_array;
//	}
//	else{
//		printf("Unbkown vector field\n");
//		return 0;
//	}
//}


template<int dimension>
Array2D<2> * Array_ZoomIn<dimension>::generateVectorField(VectorFieldMode vector_field)
{
	Array2D<2> * vector_field_array = new Array2D<2>(height, width, DUAL_DUAL, MIRROR_EXTENSION);
#pragma omp parallel for
		for (int j = 0; j < height; j++){
			for (int i = 0; i < width; i++){
				Point<2> sample_position = vector_field_array->arrayIndexTextureCoordinate(i, j);
				if (vector_field == MEAN_GRADIENT_FIELD){
					sample_position = currentFieldToInitalFieldCoordinates(sample_position);
					vector_field_array->entries[i + j * width] = MeanGradient(sample_position);
				}
				else if (vector_field == ADVECTION_VECTOR_FIELD){
					sample_position = currentFieldToInitalFieldCoordinates(sample_position);
					vector_field_array->entries[i + j * width] = AdvectionVector(sample_position);
				}
				else if (vector_field == ACCUMULATION_VECTOR_FIELD){
					vector_field_array->entries[i + j * width] = Point<2>(gradient_accumulation_w->applyContinousFilter(sample_position).Mean(), gradient_accumulation_h->applyContinousFilter(sample_position).Mean());
				}
				//else if (vector_field == HESSIAN_PRINCIPAL_FIELD){
				//	Matrix2D hessian = MeanHessianMatrix(sample_position);
				//	Point<2> principal_directions[2];
				//	double eigenvalues[2];
				//	Matrix2DPrincipalComponents(hessian, principal_directions, eigenvalues);
				//	vector_field_array->entries[i + j * width] = principal_directions[0] * eigenvalues[0];
				//}
				//else if (vector_field == HESSIAN_SECONDARY_FIELD){
				//	Matrix2D hessian = MeanHessianMatrix(sample_position);
				//	Point<2> principal_directions[2];
				//	double eigenvalues[2];
				//	Matrix2DPrincipalComponents(hessian, principal_directions, eigenvalues);
				//	vector_field_array->entries[i + j * width] = principal_directions[1] * eigenvalues[1];
				//}
				else{
					printf("Unkown vector field mode \n");
				}
			}
		}
		return vector_field_array;
}


template<int dimension>
void Array_ZoomIn<dimension>::Advect()
{
	double start_time = GetTime();
	Array2D<2> * current_end_position = PathAdvection();
	//Array2D<2> * current_end_position = ProgressiveAdvection();
	printf("Advection Time = %.4f(s)\n", GetTime() - start_time);

	UpdateGradientAccumulation(current_end_position);
	delete current_end_position;

//	printf("ADVECTION \n");
//	Array2D<2> * advection_field = generateVectorField(ADVECTION_VECTOR_FIELD);
//	advection_field->statistics->verbose = false;
//	advection_field->statistics->setArrayStatisticsL2();
//
//	double advection_field_scale = advection_field->statistics->max_value_l2_norm;
//	delete advection_field;
//
//
//	Array2D<dimension> * accumulation_w;
//	Array2D<dimension> * accumulation_h;
//	printf("GRADIENT ADVECTION \n");
//	accumulation_w = new Array2D<dimension>(height, width + 1, PRIMAL_DUAL, extension_mode);
//	accumulation_h = new Array2D<dimension>(height + 1, width, DUAL_PRIMAL, extension_mode);
//
//	bool ortogonal_path;
//	ortogonal_path = path_mode == ORTHOGONAL_PATH ? true : false;
//
//#pragma omp parallel for
//	for (int j = 0; j < height; j++){
//		for (int i = 0; i < width; i++){
//		
//			Point<2> current_field_initial_position = current_field->arrayIndexTextureCoordinate(i, j);
//			Point<2> initial_field_initial_position = currentFieldToInitalFieldCoordinates(current_field_initial_position);
//			Point<2> advection_direction = AdvectionVector(initial_field_initial_position);
//
//			//Point<dimension> gradient_w = current_field->applyDerivativeFilter(initial_position, DX);
//			//Point<dimension> gradient_h = current_field->applyDerivativeFilter(initial_position, DY);
//
//			Point<dimension> gradient_w = gradient_accumulation_w->applyContinousFilter(current_field_initial_position);
//			Point<dimension> gradient_h = gradient_accumulation_h->applyContinousFilter(current_field_initial_position);
//
//			double gradient_norm = gradient_w.norm() + gradient_h.norm();
//
//			if (advection_direction.norm() > 0.f && gradient_norm > 0.f){ // Invariance:  normalized_last_advection_direction != zero-vector
//				Point<2> normalized_advection_direction = advection_direction.normalvector();
//				PathData initial_data(initial_field_initial_position, advection_time / advection_field_scale, normalized_advection_direction);
//
//				PathData forward_data;
//				forward_data = PathIntegration(initial_data, false, ortogonal_path);
//
//#pragma omp critical
//				{
//					TransportFunction_Gradient(gradient_w, gradient_h, normalized_advection_direction,
//					initalFieldToCurrentFieldCoordinates(forward_data.current_position), forward_data.normalized_last_advection_direction,
//					1.f, 1.f, 1.f, 1.f, 1, accumulation_w, accumulation_h);
//				}
//			}
//			else if(gradient_norm > 0.f){
//#pragma omp critical
//				{
//					accumulation_w->splatValue(current_field_initial_position, gradient_w);
//					accumulation_h->splatValue(current_field_initial_position, gradient_h);
//				}
//			}
//		}
//	}
//
//	gradient_accumulation_w->set(accumulation_w);
//	gradient_accumulation_h->set(accumulation_h);
//
//	delete accumulation_w;
//	delete accumulation_h;
//
//	printf("ADVECTION DONE!!\n");

}

template<int dimension>
void Array_ZoomIn<dimension>::TransportFunction_Gradient(Point<dimension> gradient_w, Point<dimension> gradient_h, Point<2> normalized_initial_advection_direction,
	Point<2> current_position, Point<2> normalized_last_advection_direction,
	double time_difference, double spatial_difference, double path_total_time, double path_total_lenght, int num_divisions,
	Array2D<dimension> * accumulation_w, Array2D<dimension> * accumulation_h)
{
	double weight;
	if (num_divisions == 1){
		weight = 1.f;
	}
	else{
		if (accumulation_weighting_mode == DIRECT_SPATIAL_WEIGHTING){
			weight = spatial_difference / path_total_lenght;
		}
		else if (accumulation_weighting_mode == INVERSE_SPATIAL_WEIGHTING){
			weight = (path_total_lenght - spatial_difference) / (static_cast<double>(num_divisions - 1)*path_total_lenght);
		}
		else if (accumulation_weighting_mode == DIRECT_TEMPORAL_WEIGHTING){
			weight = time_difference / path_total_time;
		}
		else if (accumulation_weighting_mode == INVERSE_TEMPORAL_WEIGHTING){
			weight = (path_total_time - time_difference) / (static_cast<double>(num_divisions - 1)*path_total_time);
		}
	}

	if (transport_mode == PARALLEL_TRANSPORT){
		Point<dimension> rotated_gradient_w;
		Point<dimension> rotated_gradient_h;

		//Point<dimension> current_gradient_w = current_field->applyDerivativeFilter(current_position,DX);
		//Point<dimension> current_gradient_h = current_field->applyDerivativeFilter(current_position,DY);

		//Point<dimension> current_gradient_w = gradient_accumulation_w->applyContinousFilter(current_position);
		//Point<dimension> current_gradient_h = gradient_accumulation_h->applyContinousFilter(current_position);

		vectorPlanarRotation(gradient_w, gradient_h, rotated_gradient_w, rotated_gradient_h, normalized_initial_advection_direction, normalized_last_advection_direction);
		accumulation_w->splatValue(current_position, rotated_gradient_w*weight);
		accumulation_h->splatValue(current_position, rotated_gradient_h*weight);
	}
	else if (transport_mode == MAGNITUDE_TRANSPORT)
	{
		//Point<dimension> current_gradient_w = current_field->applyDerivativeFilter(current_position,DX);
		//Point<dimension> current_gradient_h = current_field->applyDerivativeFilter(current_position,DY);

		Point<dimension> current_gradient_w = gradient_accumulation_w->applyContinousFilter(current_position);
		Point<dimension> current_gradient_h = gradient_accumulation_h->applyContinousFilter(current_position);

		double current_gradient_norm = sqrt(current_gradient_w.squared_norm() + current_gradient_h.squared_norm());
		if (current_gradient_norm > 0.f){
			double initial_gradient_norm = sqrt(gradient_w.squared_norm() + gradient_h.squared_norm());
			double scale_factor = initial_gradient_norm / current_gradient_norm;
			accumulation_w->splatValue(current_position, current_gradient_w*scale_factor*weight);
			accumulation_h->splatValue(current_position, current_gradient_h*scale_factor*weight);
		}
		else{
			accumulation_w->splatValue(current_position, gradient_w*weight);
			accumulation_h->splatValue(current_position, gradient_h*weight);
		}
	}
	else if (transport_mode == DIRECT_TRANSPORT){
		accumulation_w->splatValue(current_position, gradient_w*weight);
		accumulation_h->splatValue(current_position, gradient_h*weight);
	}
}

template<int dimension>
PathData  Array_ZoomIn<dimension>::PathIntegration(PathData input, bool backward_direction, bool orthogonal_direction)
{
	int function_calls = 0;
	while (function_calls<max_function_calls && input.time > 0.f)
	{
		Point<2> advection_direction = AdvectionVector(input.current_position);
		if (orthogonal_direction){
			advection_direction = Point<2>(-advection_direction[1], advection_direction[0]);
		}

		double advection_direction_norm = advection_direction.norm();
		if (advection_direction_norm == 0.f){
			return input;
		}
		Point<2> normalized_advection_direction = advection_direction / advection_direction_norm;
		if (normalized_advection_direction.dot(input.normalized_last_advection_direction) < max_angle_deviation){// Test with coplanar directions
			return input;
		}
		double advection_time = input.time;
		if (advection_direction_norm*input.time > max_step_lenght){
			advection_time = (max_step_lenght / advection_direction_norm);
		}

		Point<2> new_position = backward_direction ? input.current_position - advection_direction*advection_time : input.current_position + advection_direction*advection_time;
		if (new_position[0] < 0.f || new_position[0] > static_cast<double>(width) || new_position[1] < 0.f || new_position[1] > static_cast<double>(height)){
			return input;
		}

		input.current_position = new_position;
		input.time -= advection_time;
		input.normalized_last_advection_direction = normalized_advection_direction;
		input.lenght += advection_direction_norm*advection_time;
		input.divisions++;

		function_calls++;
	}
	if (function_calls == max_function_calls){
		printf("Maximum of function calls attained !!!\n");
	}
	return input;
}

//template<int dimension>
//void Array_ZoomIn<dimension>::UpdateGradientAccumulation()
//{
//	printf("UPDATING GRADIENT \n");
//	gradient_accumulation_w->getValuesFromDiscreteFiltering(current_field, derivative_w);
//	gradient_accumulation_w->statistics->updateStatistics();
//	gradient_accumulation_h->getValuesFromDiscreteFiltering(current_field, derivative_h);
//	gradient_accumulation_h->statistics->updateStatistics();
//	printf("DONE!! \n");
//}

//template<int dimension>
//void Array_ZoomIn<dimension>::UpdateMeanCurrentField()
//{
//	printf("UPDATING MEAN CURRENT FIELD \n");
//	current_field->MeanArray(mean_initial_field);
//	printf("DONE!! \n");
//}

//template<int dimension>
//void Array_ZoomIn<dimension>::UpdateCurrentFieldDependants()
//{
//	UpdateGradientAccumulation();
//	UpdateMeanCurrentField();
//}


template<int dimension>
SparseSystem<dimension> * Array_ZoomIn<dimension>::SetLinearSystem(LinearSystemMode system_mode)
{
	int num_discrete_operators;
	DiscreteFilter ** discrete_operators;
	Array2D<dimension> ** operator_values;
	double * weights;
	double * scaling;

	if (system_mode == GRADIENT_FIELD_LEAST_SQUARES)
	{
		num_discrete_operators = 2;
		discrete_operators = new  DiscreteFilter *[num_discrete_operators];
		operator_values = new  Array2D<dimension> *[num_discrete_operators];
		weights = new double[num_discrete_operators];
		scaling = new double[num_discrete_operators];

		discrete_operators[0] = derivative_w;
		operator_values[0] = gradient_accumulation_w;

		discrete_operators[1] = derivative_h;
		operator_values[1] = gradient_accumulation_h;

		weights[0] = weights[1] = 1.f;
		scaling[0] = scaling[1] = 1.f;
	}

	SparseSystem<dimension> * linear_system = SetLinearSystemCoefficients(num_discrete_operators, discrete_operators, operator_values, weights, scaling);
	delete discrete_operators;
	delete operator_values;
	delete weights;

	printf("DONE \n");
	if (linear_system->CheckSystemSymmetry()){
		printf("SYMMETRIC SYSTEM\n");
	}
	else{
		printf("NON SYMMETRIC SYSTEM!! \n");
	}
	if (linear_system->CheckDiagonalDominance()){
		printf("DIAGONAL DOMINANT SYSTEM\n");
	}
	else{
		printf("NON DIAGONAL DOMINANT SYSTEM!! \n");
	}

	return linear_system;
}

template<int dimension>
SparseSystem<dimension> * Array_ZoomIn<dimension>::SetLinearSystemCoefficients(int num_discrete_operators, DiscreteFilter ** discrete_operators, Array2D<dimension> ** operator_values, double * weights, double * scaling)
{
	SparseSystem<dimension> * linear_system = new SparseSystem<dimension>(current_field->array_size);

	int system_radius_h = 0;
	int system_radius_w = 0;
	for (int iter = 0; iter < num_discrete_operators; iter++){
		system_radius_h = discrete_operators[iter]->radius_h > system_radius_h ? discrete_operators[iter]->radius_h : system_radius_h;
		system_radius_w = discrete_operators[iter]->radius_w > system_radius_w ? discrete_operators[iter]->radius_w : system_radius_w;
	}

	double temp[dimension];
	// Traverse Array
//#pragma omp parallel for
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			// Traverse Neighbourhood
			for (int dj = -system_radius_h; dj < system_radius_h + 1; dj++){
				for (int di = -system_radius_w; di < system_radius_w + 1; di++){
					int pi = ExtensionIndex(i + di, width, extension_mode);
					int pj = ExtensionIndex(j + dj, height, extension_mode);
					int current_index = pi + pj*width;

					for (int cj = -system_radius_h; cj < system_radius_h + 1; cj++){
						for (int ci = -system_radius_w; ci < system_radius_w + 1; ci++){

							int qi = ExtensionIndex(i + ci, width, extension_mode);
							int qj = ExtensionIndex(j + cj, height, extension_mode);
							int neighbour_index = qi + qj*width;

							double coefficient = 0.f;
							for (int operator_index = 0; operator_index < num_discrete_operators; operator_index++){
								coefficient += discrete_operators[operator_index]->value(di, dj)*discrete_operators[operator_index]->value(ci, cj)*weights[operator_index];
							}
//#pragma omp critical
							{
								linear_system->elements[current_index].updateCoefficient(neighbour_index, coefficient);
							}
							
						}
					}

					Point<dimension> constraint;
					for (int operator_index = 0; operator_index < num_discrete_operators; operator_index++){
						constraint += (operator_values[operator_index]->value(i, j)  * scaling[operator_index])* discrete_operators[operator_index]->value(di, dj) * weights[operator_index];
					}

					for (int k = 0; k < dimension; k++){ temp[k] = constraint[k]; }
//#pragma omp critical
					{
						linear_system->elements[current_index].updateConstraint(temp);
					}
					
				}
			}
			for (int k = 0; k < dimension; k++){
				linear_system->elements[i + j*width].initial_value[k] = current_field->value(i, j)[k];
			}

		}
	}
	return linear_system;
}

template<int dimension>
void Array_ZoomIn<dimension>::SolveLinearSystem(LinearSystemMode system_mode, bool correct_dc_solution, bool clamp_solution)
{
	printf("SETTING LINEAR SYSTEM \n");
	SparseSystem<dimension> * linear_system = SetLinearSystem(system_mode);
	printf("DONE \n");

	printf("SOLVING LINEAR SYSTEM \n");
	int cg_iterations = 100;
	printf("Max CG Iterations = %d \n", cg_iterations);
	double epsilon = 0.000001f;
	printf("CG Epsilon = %g \n", epsilon);

	//bool correct_dc_solution;
	//bool clamp_solution;

	//setSystemSolutionParameters(correct_dc_solution, clamp_solution);

	double * solution = linear_system->CGSolver(epsilon, cg_iterations, correct_dc_solution, clamp_solution);
	printf("DONE \n");
	current_field->getValuesFromExternal(solution);
	delete solution;
	delete linear_system;

	
	//UpdateCurrentFieldDependants();

}

template<int dimension>
Point<2>  Array_ZoomIn<dimension>::currentFieldToInitalFieldCoordinates(Point<2> current_coordinates)
{
	current_coordinates *= current_scale;
	current_coordinates += current_translation;

	return current_coordinates;
}

template<int dimension>
Point<2>  Array_ZoomIn<dimension>::initalFieldToCurrentFieldCoordinates(Point<2> initial_coordinates)
{
	initial_coordinates -= current_translation;
	initial_coordinates /= current_scale;

	return initial_coordinates;
}

template<int dimension>
void Array_ZoomIn<dimension>::ExtractCurrentState(){
	printf("EXTRACTION \n");

	//int corner_w = static_cast<int>(static_cast<double>(width)*zoom_box->canvas_coordinate[0]);
	//int corner_h = static_cast<int>(static_cast<double>(height)*zoom_box->canvas_coordinate[1]);

	//int cropping_width = static_cast<int>(floor(static_cast<double>(current_field->width) / zoom_box->zoom_factor));
	//int cropping_height = static_cast<int>(floor(static_cast<double>(current_field->height) / zoom_box->zoom_factor));

	//current_translation += Point<>

	//Array2D<dimension> * cropped_current_field = cropArray(current_field, corner_w, corner_h, cropping_width, cropping_height);
	//current_field->getValuesFromContinuousFiltering(cropped_current_field,NORMALIZED_TEXTURE_COORD);

	//delete cropped_current_field;

	double corner_w = static_cast<double>(width)*zoom_box->canvas_coordinate[0];
	double corner_h = static_cast<double>(height)*zoom_box->canvas_coordinate[1];

	current_translation += Point<2>(corner_w, corner_h)*current_scale;

	double region_width = static_cast<double>(width) / zoom_box->zoom_factor;
	double region_height = static_cast<double>(height) / zoom_box->zoom_factor;

	current_scale /= zoom_box->zoom_factor;

	current_field->RegionUpsampling(corner_w, corner_h, region_width, region_height, 1.f, zoom_filter);

	gradient_accumulation_w->RegionUpsampling(corner_w, corner_h, region_width, region_height, 1.f / zoom_box->zoom_factor, zoom_filter);
	gradient_accumulation_h->RegionUpsampling(corner_w, corner_h, region_width, region_height, 1.f / zoom_box->zoom_factor, zoom_filter);

	//SampleGradientAccumulationFromCurrentField();

	max_step_lenght /= zoom_box->zoom_factor;

	//Array2D<dimension> * cropped_current_accumulation_w = cropArray(current_accumulation_w, corner_w, corner_w, cropping_width, cropping_height);
	//current_accumulation_w->getValuesFromContinuousFiltering(cropped_current_accumulation_w);
	//delete cropped_current_accumulation_w;

	//Array2D<dimension> * cropped_current_accumulation_h = cropArray(current_accumulation_h, corner_w, corner_w, cropping_width, cropping_height);
	//current_accumulation_h->getValuesFromContinuousFiltering(cropped_current_accumulation_h);
	//delete cropped_current_accumulation_h;

	Point<3> canvas_coordinate = Point<3>(0.25f, 0.25f, (canvas->edges_lenght[0] + canvas->edges_lenght[1]) / 1000.f);
	zoom_box->canvas_coordinate = canvas_coordinate;

	printf("DONE!! \n");
}


template<int dimension>
void Array_ZoomIn<dimension>::SampleZoomRegion(){

	double corner_w = zoom_box->canvas_coordinate[0];
	double corner_h = zoom_box->canvas_coordinate[1];

	double region_width = 1.f / zoom_box->zoom_factor;
	double region_height = 1.f / zoom_box->zoom_factor;

	double start_time = GetTime();
	//current_field->getValuesFromRegionUpsampling(initial_field,corner_w, corner_h, region_width, region_height, 1.f, zoom_filter);
	current_field->getValuesFromRegionUpsampling(initial_field, corner_w, corner_h, region_width, region_height, 1.f);
	//printf("Sampling Time %.4f(s) \n", GetTime() - start_time);

	pixel_texture_buffer_updated = false;
}


template<int dimension>
void Array_ZoomIn<dimension>::Advance()
{
	ExtractCurrentState();
	Advect();
	SolveLinearSystem(GRADIENT_FIELD_LEAST_SQUARES);
}

template<int dimension>
Array2D<2> * Array_ZoomIn<dimension>::ProgressiveAdvection(int iterations, int clamping_percentile, double field_scale, ContinuousFilter * filter)
{
	Array2D<2> * vector_field = generateVectorField(ADVECTION_VECTOR_FIELD);
	vector_field->statistics->verbose = true;
	vector_field->statistics->setArrayStatisticsL2();

	double clamping_value = vector_field->statistics->value_l2_norm_percentiles[clamping_percentile];
	printf("Clamping Value %g \n", clamping_value);
	Array2D<2> * position_field = new Array2D<2>(height, width, DUAL_DUAL, MIRROR_EXTENSION);

#pragma omp parallel for
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			double vector_field_norm = vector_field->value(i, j).norm();
			double scalar_factor = vector_field_norm > clamping_value ? vector_field_norm : clamping_value;
			position_field->value(i, j) = position_field->arrayIndexTextureCoordinate(i, j) + vector_field->value(i, j)*field_scale / scalar_factor;
		}
	}

	delete vector_field;
	for (int d = 0; d < iterations; d++){
		UpdatePositionField(position_field, filter);
	}

	return position_field;
}

template<int dimension>
void Array_ZoomIn<dimension>::UpdatePositionField(Array2D<2> * position_field, ContinuousFilter * filter)
{

	Array2D<2> * new_position_field = new Array2D<2>(height, width, DUAL_DUAL, MIRROR_EXTENSION);
#pragma omp parallel for
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			new_position_field->value(i, j) = position_field->applyContinousFilter(position_field->value(i, j), filter);
		}
	}
	position_field->set(new_position_field);
	delete new_position_field;
}

template<int dimension>
Array2D<2> * Array_ZoomIn<dimension>::PathAdvection()
{
	printf("PATH ADVECTION \n");
	Array2D<2> * advection_field = generateVectorField(ADVECTION_VECTOR_FIELD);
	advection_field->statistics->verbose = false;
	advection_field->statistics->setArrayStatisticsL2();

	double advection_field_scale = advection_field->statistics->max_value_l2_norm;
	delete advection_field;

	bool ortogonal_path;
	ortogonal_path = path_mode == ORTHOGONAL_PATH ? true : false;

	Array2D<2> * current_end_position = new Array2D<2>(height, width, DUAL_DUAL, MIRROR_EXTENSION);

#pragma omp parallel for
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){

			Point<2> current_field_initial_position = current_field->arrayIndexTextureCoordinate(i, j);
			Point<2> initial_field_initial_position = currentFieldToInitalFieldCoordinates(current_field_initial_position);
			Point<2> advection_direction = AdvectionVector(initial_field_initial_position);

			if (advection_direction.norm() > 0.f){ // Invariance:  normalized_last_advection_direction != zero-vector
				Point<2> normalized_advection_direction = advection_direction.normalvector();
				PathData initial_data(initial_field_initial_position, advection_time / advection_field_scale, normalized_advection_direction);

				PathData forward_data;
				forward_data = PathIntegration(initial_data, false, ortogonal_path);
				current_end_position->value(i, j) = initalFieldToCurrentFieldCoordinates(forward_data.current_position);
			}
		}
	}
	return current_end_position;
	printf("DONE!!\n");
}

template<int dimension>
void Array_ZoomIn<dimension>::UpdateGradientAccumulation(Array2D<2> * current_advected_position)
{
	Array2D<dimension> * accumulation_w;
	Array2D<dimension> * accumulation_h;

	printf("UPDATE ACCUMULATION \n");
	accumulation_w = new Array2D<dimension>(height, width + 1, PRIMAL_DUAL, extension_mode);
	accumulation_h = new Array2D<dimension>(height + 1, width, DUAL_PRIMAL, extension_mode);

#pragma omp parallel for
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			Point<2> current_start_position = current_field->arrayIndexTextureCoordinate(i, j);
			Point<2> current_end_position = current_advected_position->value(i, j);
#pragma omp critical
				{
					GradientAccumulation(current_start_position,current_end_position,accumulation_w,accumulation_h);
				}
		}
	}

	gradient_accumulation_w->set(accumulation_w);
	gradient_accumulation_h->set(accumulation_h);

	delete accumulation_w;
	delete accumulation_h;

	printf("DONE!!\n");

}
template<int dimension>
void Array_ZoomIn<dimension>::GradientAccumulation(Point<2> current_start_position, Point<2> current_end_position, Array2D<dimension> * accumulation_w, Array2D<dimension> * accumulation_h)
{
	Point<dimension> gradient_w = gradient_accumulation_w->applyContinousFilter(current_start_position);
	Point<dimension> gradient_h = gradient_accumulation_h->applyContinousFilter(current_start_position);

	double gradient_norm = gradient_w.norm() + gradient_h.norm();

	if (gradient_norm > 0.f){

		if (transport_mode == PARALLEL_TRANSPORT){
			Point<dimension> rotated_gradient_w;
			Point<dimension> rotated_gradient_h;

			Point<2> intial_start_position = currentFieldToInitalFieldCoordinates(current_start_position);
			Point<2> intial_end_position = currentFieldToInitalFieldCoordinates(current_end_position);

			Point<2> start_advection_direction = MeanGradient(intial_start_position).normalize();
			Point<2> end_advection_direction = MeanGradient(intial_end_position).normalize();

			vectorPlanarRotation(gradient_w, gradient_h, rotated_gradient_w, rotated_gradient_h, start_advection_direction, end_advection_direction);
			accumulation_w->splatValue(current_end_position, rotated_gradient_w);
			accumulation_h->splatValue(current_end_position, rotated_gradient_h);
		}
		else if (transport_mode == MAGNITUDE_TRANSPORT)
		{
			Point<dimension> current_gradient_w = gradient_accumulation_w->applyContinousFilter(current_end_position);
			Point<dimension> current_gradient_h = gradient_accumulation_h->applyContinousFilter(current_end_position);

			double current_gradient_norm = sqrt(current_gradient_w.squared_norm() + current_gradient_h.squared_norm());
			if (current_gradient_norm > 0.f){
				double initial_gradient_norm = sqrt(gradient_w.squared_norm() + gradient_h.squared_norm());
				double scale_factor = initial_gradient_norm / current_gradient_norm;
				accumulation_w->splatValue(current_end_position, current_gradient_w*scale_factor);
				accumulation_h->splatValue(current_end_position, current_gradient_h*scale_factor);
			}
			else{
				accumulation_w->splatValue(current_end_position, gradient_w);
				accumulation_h->splatValue(current_end_position, gradient_h);
			}
		}
		else if (transport_mode == DIRECT_TRANSPORT){
			accumulation_w->splatValue(current_end_position, gradient_w);
			accumulation_h->splatValue(current_end_position, gradient_h);
		}
	}
}

template<int dimension>
void Array_ZoomIn<dimension>::SetCholeskyFactorization(LinearSystemMode system_mode)
{
	SparseSystem<dimension> * linear_system = SetLinearSystem(system_mode);

	SparseMatrixd A = linear_system->GenerateSparseMatrix();
	double start_time = GetTime();
	cholesky_factorization = new CholeskyFactorization(A);
	printf("Cholesky Factorization = %.4f(s)\n", GetTime() - start_time);
	delete linear_system;
}

template<int dimension>
void Array_ZoomIn<dimension>::SolveFromCholesky(LinearSystemMode system_mode)
{
	SparseSystem<dimension> * linear_system = SetLinearSystem(system_mode);

	for (int channel = 0; channel < dimension; channel++)
	{
		printf("Channel %d \n",channel);
		Vectord rhs_vector = linear_system->GenerateRHS(channel);
		double initial_dc = linear_system->InitialDC(channel);
		printf("Initial DC %f \n", initial_dc);
		double start_time = GetTime();
		Vectord solution = cholesky_factorization->solve(rhs_vector);
		printf("Cholesky Solution = %.4f(s)\n", GetTime() - start_time);
		double final_dc = solution.sum();
		printf("Final DC %f \n", final_dc);
		current_field->getValuesFromVectord(solution, channel, ((initial_dc - final_dc) / static_cast<double>(array_size)));
	}

	delete linear_system;
}