#include "geometry.h"

Matrix3D ApplyNeighbourTransformation(Matrix3D transformation_to_neighbour, Matrix3D input, Point<3> normal, bool apply_transformation)
{
	if (apply_transformation)
	{
		Matrix3D output;
		for (int i = 0; i < 3; i++)
		{
			Point<3> channel_gradient = input.ExtractRow(i);
			Point<3> transformed_gradient = transformation_to_neighbour * channel_gradient;
			//Project To plane to avoid error accumuluation
			Point<3> projected_gradient = transformed_gradient - normal*transformed_gradient.dot(normal);
			output.InsertRow(projected_gradient, i);
		}
		return output;
	}
	else
		return input;
}

Point<3> ApplyNeighbourTransformation(Matrix3D transformation_to_neighbour, Point<3> input, Point<3> normal, bool apply_transformation)
{
	if (apply_transformation)
	{
		return transformation_to_neighbour * input;
	}
	else
		return input;

}

Point<3> MinPrincipalDirection(Matrix3D correlation, Point<3> x0 = Point<3>(1.f,0.f,0.f))
{
	double epsilon = 0.00001f;
	double max_step = 0.1f;
	int max_iter = 20;

	Point<3> x_new = x0;
	Point<3> x_prev;
	Point<3> gradient;
	Point<3> direction;
	double prev_val;
	double new_val = DBL_MAX;
	double step;
	int iter=0;
	do{
		x_prev = x_new;
		prev_val = new_val;
		gradient = correlation*x_prev;
		new_val = gradient.dot(prev_val);
		direction = (gradient - x_prev*gradient.dot(x_prev))*(-1.f);
		step = direction.dot(direction) / direction.dot(correlation*direction);
		step = step < max_step ? step : max_step;
		x_new = x_prev + direction*step;
		x_new.normalize();
		iter++;
	} while ((x_prev - x_new).norm()> epsilon && iter < max_iter);

	return x_new;
}