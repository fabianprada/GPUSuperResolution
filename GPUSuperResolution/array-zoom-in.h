#ifndef ARRAY_ZOOMIN_INCLUDE
#define ARRAY_ZOOMIN_INCLUDE

#include "array-differential.h"
#include "zoom.h"
#include "scene.h"

template<int dimension>
class Array_ZoomIn : public GraphicObject{
public:
	Array_ZoomIn(Array2D<dimension> * p_initial_field, Camera * p_camera){

		sprintf(object_name, "ARRAY ZOOM IN");
		is_active = true;

		initial_field = new Array2D<dimension>(p_initial_field);
		current_field = new Array2D<dimension>(initial_field);

		height = initial_field->height;
		width = initial_field->width;
		array_size = initial_field->array_size;
		extension_mode  = initial_field->array_extension_mode;

		mean_initial_field = new Array2D<1>(height,width, initial_field->sampling_mode, extension_mode);
		current_field->MeanArray(mean_initial_field);
		//gradient_field_w = new Array2D<dimension>(height, width + 1, PRIMAL_DUAL, extension_mode);
		//gradient_field_h = new Array2D<dimension>(height + 1, width, DUAL_PRIMAL, extension_mode);

		gradient_accumulation_w = new Array2D<dimension>(height, width + 1, PRIMAL_DUAL, extension_mode);
		gradient_accumulation_h = new Array2D<dimension>(height + 1, width, DUAL_PRIMAL, extension_mode);

		derivative_h = new Derivative_h();
		derivative_w = new Derivative_w();
		SampleGradientAccumulationFromCurrentField();

		//UpdateCurrentFieldDependants();
		
		initial_texture_canvas = new Canvas_Rectangle(Point<3>(-1.f, -static_cast<double>(height) / static_cast<double>(width), 0.f), Point<3>(1.f, 0.f, 0.f),
		Point<3>(0.f, 1.f, 0.f), 2.f, 2.f*static_cast<double>(height) / static_cast<double>(width));
		initial_texture_canvas->reference_vertex += Point<3>(-1.2f, 0.f, 0.f);
		Point<3> canvas_coordinate = Point<3>(0.25f, 0.25f, (initial_texture_canvas->edges_lenght[0] + initial_texture_canvas->edges_lenght[1]) / 1000.f);
		zoom_box = new ZoomBox(initial_texture_canvas, canvas_coordinate);


		canvas = new Canvas_Rectangle(Point<3>(-1.f, -static_cast<double>(height) / static_cast<double>(width), 0.f), Point<3>(1.f, 0.f, 0.f),
		Point<3>(0.f, 1.f, 0.f), 2.f, 2.f*static_cast<double>(height) / static_cast<double>(width));
		canvas->reference_vertex += Point<3>(1.2f, 0.f, 0.f);
		//Point<3> canvas_coordinate = Point<3>(0.25f, 0.25f, (canvas->edges_lenght[0] + canvas->edges_lenght[1]) / 1000.f);
		//zoom_box = new ZoomBox(canvas, canvas_coordinate);

		camera = p_camera;

		zoom_filter = new ContinuousBspline3();

		vector_visualization_scale_factor = 2.f;

		visualization_mode = CURRENT_VALUE_VISUALIZATION;
		vector_field_mode = NONE_VECTOR_FIELD;
		advection_vector_mode = NORMALIZED_MEAN_GRADIENT_PROJECTED_GRADIENT;
		transport_mode = DIRECT_TRANSPORT;
		accumulation_weighting_mode = INVERSE_SPATIAL_WEIGHTING;
		path_mode = PARALLEL_PATH;

		advection_time = 100.f;
		max_angle_deviation = 0.1f;
		max_step_lenght = 0.5f;
		max_function_calls = 200;

		current_translation = Point<2>();
		current_scale = 1.f;

		SampleZoomRegion();

		//SetCholeskyFactorization(GRADIENT_FIELD_LEAST_SQUARES);

	}

	Canvas_Rectangle * canvas;
	Canvas_Rectangle * initial_texture_canvas;

	ZoomBox * zoom_box;
	Camera * camera;

	Point<2> current_field_corner;
	double current_field_heigth;
	double current_field_width;

	Array2D<dimension>* initial_field;
	Array2D<dimension>* current_field;
	Array2D<1> * mean_initial_field;

	Point<2> current_translation;
	double current_scale;

	//Array2D<dimension> * gradient_field_w;
	//Array2D<dimension> * gradient_field_h;

	DiscreteFilter * derivative_h;
	DiscreteFilter * derivative_w;

	//void UpdateGradientAccumulation();
	//void UpdateMeanCurrentField();
	//void UpdateCurrentFieldDependants();
	void SampleGradientAccumulationFromCurrentField();
	
	Array2D<dimension>* gradient_accumulation_w;
	Array2D<dimension>* gradient_accumulation_h;

	ContinuousFilter * zoom_filter;

	void RestartGradientAccumulation();
	void SampleZoomRegion();
	void ExtractCurrentState();
	void ExtractCurrentState_CallBack(Prompt  & prompt);
	//void AdvectCurrentState();
	//void IntegrateCurrentState();

	GLuint initial_texture_buffer;
	void setupInitialTexture(GLuint & texture_buffer);

	VisualizationArray visualization_mode;
	GLuint pixel_texture_buffer;
	bool pixel_texture_buffer_updated;
	void setupTexture(GLuint & texture_buffer);
	void setTextureValues(GLuint & texture_buffer);
	void ChangeVisualizationArray_CallBack(Prompt & prompt);
	//void ChangeVisualizationArrayMisha_CallBack(Prompt & prompt);

	VectorFieldMode vector_field_mode;
	GLuint vector_field_buffer;
	bool vector_field_buffer_updated;
	void setupVectorField(GLuint & vector_buffer);
	void setVectorFieldValues(GLuint & vector_buffer);
	void ChangeVectorField_CallBack(Prompt & prompt);
	double vector_visualization_scale_factor;

	void setupOpenGL();
	void drawOpenGL(Scene * scene = 0);
	void HandleKeyboardEvent(Prompt  & prompt, Scene * scene = 0);
	void HandleMouseEvent(MouseEvent  & mouse_event, Scene * scene = 0);

	Point<2> MeanGradient(Point<2> sample_position);
	Matrix2D MeanHessian(Point<2> sample_position);
	Point<2> AdvectionVector(Point<2> sample_position);
	Array2D<2> * generateVectorField(VectorFieldMode vector_field);
	//Array2D<2> * generateDisplacementField(VectorFieldMode vector_field)

	int height, width, array_size;
	ExtensionMode extension_mode;

	double advection_time;
	double max_angle_deviation;
	double max_step_lenght;
	int max_function_calls;

	AccumulationWeightingMode accumulation_weighting_mode;
	TransportMode transport_mode;
	PathMode path_mode;
	void Advection_CallBack(Prompt & prompt);
	AdvectionVectorMode advection_vector_mode;

	PathData PathIntegration(PathData input, bool backward_direction = false, bool orthogonal_direction = false);

	Point<2> currentFieldToInitalFieldCoordinates(Point<2> current_coordinates);
	Point<2> initalFieldToCurrentFieldCoordinates(Point<2> initial_coordinates);
	
	Array2D<2> * ProgressiveAdvection(int iterations = 6, int clamping_percentile = 8, double field_scale = 0.25f, ContinuousFilter * filter = 0);
	Array2D<2> * PathAdvection();
	void UpdateGradientAccumulation(Array2D<2> * destination_position);
	void UpdatePositionField(Array2D<2> * position_field, ContinuousFilter * filter);
	void Advance_CallBack(Prompt & prompt);
	void Advance();
	void Advect();

	void GradientAccumulation(Point<2> current_start_position, Point<2> current_end_position,Array2D<dimension> * accumulation_w, Array2D<dimension> * accumulation_h);

	void TransportFunction_Gradient(Point<dimension> gradient_w, Point<dimension> gradient_h, Point<2> normalized_initial_advection_direction,
		Point<2> current_position, Point<2> normalized_last_advection_direction,
		double time_difference, double spatial_difference, double path_total_time, double path_total_lenght, int num_divisions,
		Array2D<dimension> * accumulation_w, Array2D<dimension> * accumulation_h);

	//Linear System
	void IntegrateGradientField_CallBack(Prompt & prompt);
	void SolveLinearSystem(LinearSystemMode system_mode, bool correct_dc_solution = true, bool clamp_solution = true);
	SparseSystem<dimension> * SetLinearSystem(LinearSystemMode system_mode);
	SparseSystem<dimension> * SetLinearSystemCoefficients(int num_discrete_operators, DiscreteFilter ** discrete_operators, Array2D<dimension> ** operator_values, double * weights, double * scaling);

	void SetCholeskyFactorization(LinearSystemMode system_mode);
	CholeskyFactorization * cholesky_factorization;
	void SolveFromCholesky(LinearSystemMode system_mode);
};

#include "array-zoom-in-rendering.inl"
#include "array-zoom-in.inl"
#endif //ARRAY_ZOOMIN_INCLUDE