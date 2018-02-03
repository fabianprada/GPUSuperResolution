#include "geometry.h"
#include "blinn.h"

void Matrix2DPrincipalComponents(Matrix2D matrix, Point<2> * principal_directions, double * eigenvalues)
{
	if (matrix.m[1][0] == matrix.m[0][1]){
		double a = 1.f;
		double b = -(matrix.m[0][0] + matrix.m[1][1]);
		double c = (matrix.m[0][0] * matrix.m[1][1] - matrix.m[0][1] * matrix.m[1][0]);
		double * r = new double[2];

		if (quadraticsolve(a, b, c, r) == 0){
			printf("Non real solutions!!\n");
		}
		else{
			double temp;
			if (abs(r[0]) < abs(r[1])){ temp = r[0]; r[0] = r[1]; r[1] = temp; }
			eigenvalues[0] = r[0];
			eigenvalues[1] = r[1];
			if (r[0] != r[1]){
				double * solution = new double[2];
				if (solve2x2(matrix.m[0][0] - r[0], matrix.m[0][1], matrix.m[1][0], matrix.m[1][1] - r[0], solution) != 0){
					principal_directions[0] = Point<2>(solution[0], solution[1]);
					principal_directions[1] = Point<2>(-solution[1], solution[0]);
				}
			}
			else
			{
				printf(" Umbilical Point r0 =%g,  r1 =%g \n", r[0], r[1]);
				principal_directions[0] = Point<2>(1.f, 0.f);
				principal_directions[1] = Point<2>(0.f, 1.f);
			}
		}
	}
	else{
		printf("Non symmetric matrix!!\n");
	}
}

void Matrix2DPrincipalComponents(Matrix2D matrix, Point<2> & principal_direction, double & principal_eigenvalue, double & secondary_eigenvalue)
{
	principal_direction.coord[0] = 1.f;
	principal_direction.coord[1] = 0.f;
	principal_eigenvalue = 0.f;
	secondary_eigenvalue = 0.f;

	if (matrix.m[1][0] == matrix.m[0][1]){
		double a = 1.f;
		double b = -(matrix.m[0][0] + matrix.m[1][1]);
		double c = (matrix.m[0][0] * matrix.m[1][1] - matrix.m[0][1] * matrix.m[1][0]);
		double * r = new double[2];

		if (quadraticsolve(a, b, c, r) == 0){
			printf("Non real solutions!!\n");
		}
		else{
			double temp;
			if (abs(r[0]) < abs(r[1])){ temp = r[0]; r[0] = r[1]; r[1] = temp; }
			principal_eigenvalue = r[0];
			secondary_eigenvalue = r[1];
			if (r[0] != r[1]){
				double * solution = new double[2];
				if (solve2x2(matrix.m[0][0] - r[0], matrix.m[0][1], matrix.m[1][0], matrix.m[1][1] - r[0], solution) != 0){
					principal_direction.coord[0] = solution[0];
					principal_direction.coord[1] = solution[1];
				}
			}
			else
			{
				printf(" Umbilical Point r0 =%g,  r1 =%g \n", r[0], r[1]);
				principal_direction.coord[0] = 1.f;
				principal_direction.coord[1] = 0.f;
			}
		}
	}
	else{
		printf("Non symmetric matrix!!\n");
	}
}

//Matrix3D RankOneMatrix(Point<3> u, Point<3> v)
//{
//	Matrix3D matrix;
//	matrix.m[0][0] = u[0] * v[0]; matrix.m[1][0] = u[0] * v[1]; matrix.m[2][0] = u[0] * v[2];
//	matrix.m[0][1] = u[1] * v[0]; matrix.m[1][1] = u[1] * v[1]; matrix.m[2][1] = u[1] * v[2];
//	matrix.m[0][2] = u[2] * v[0]; matrix.m[1][2] = u[2] * v[1]; matrix.m[2][2] = u[2] * v[2];
//
//	return matrix;
//}

