#ifndef GEOMETRY_INCLUDED
#define GEOMETRY_INCLUDED

#include <math.h>
#include <stdlib.h>
#include <cfloat>
#include <vadefs.h>
#include <stdarg.h>
#include <fstream>
#include <string.h>
#include <vector>
#include <list>
#include <set>
#include <queue>

template< int dimension>
class Point
{
public:
	Point(void){ for (int i = 0; i < dimension; i++){ coord[i] = 0.f; } }
	Point(double c0, ...);
	double coord[dimension];
	double operator[](int index) const;
	//Point & operator= (const Point & pt);// This is the standard assignment operator
	Point operator+ (const Point & pt) const;
	Point & operator+= (const Point & pt);
	Point operator- (const Point & pt) const;
	Point & operator-= (const Point & pt);
	Point operator* (const double s) const;
	Point operator* (const Point & pt) const;
	Point & operator*= (const double s);
	Point & operator*= (const Point & pt);
	Point operator/ (const double s) const;
	Point operator/ (const Point & pt) const;
	Point & operator/= (const double s);
	Point & operator/= (const Point & pt);
	Point normalvector() const;
	Point & normalize();
	double norm(void) const;
	double l1_norm(void) const;
	double linf_norm(void) const;
	double squared_norm(void) const;
	double dot(const Point pt) const;
	Point cross(const Point & pt) const;
	Point projection(const Point & pt) const;
	double min_value( int & index) const;
	double min_positive_value(int & index) const;
	double max_value( int & index) const;
	void squaredRoot();
	Point<dimension> squaredRoot() const;
	double Mean() const;
};

template<int dimension >
Point<dimension >::Point(double c0, ...)
{
	coord[0] = c0;
	va_list vl;
	va_start(vl, c0);
	for (int i = 1; i<dimension; i++) coord[i] = va_arg(vl, double);
	va_end(vl);
}

template<int dimension>
double Point<dimension>::operator[](int index) const
{
	return coord[index];
}

template<int dimension>
Point<dimension> Point<dimension>::operator+(const Point<dimension> & pt) const
{
	Point<dimension> npoint;
	for (int i = 0; i < dimension; i++){ npoint.coord[i] = coord[i] + pt.coord[i]; }
	return npoint;
}

template<int dimension>
Point<dimension>& Point<dimension>::operator+=(const Point<dimension> & pt)
{
	for (int i = 0; i < dimension; i++){ coord[i] += pt.coord[i]; }
	return *this;
}

template<int dimension>
Point<dimension> Point<dimension>::operator-(const Point<dimension> & pt) const
{
	Point<dimension> npoint;
	for (int i = 0; i < dimension; i++){ npoint.coord[i] = coord[i] - pt.coord[i]; }
	return npoint;
}

template<int dimension>
Point<dimension>& Point<dimension>::operator-=(const Point<dimension> & pt)
{
	for (int i = 0; i < dimension; i++){ coord[i] -= pt.coord[i]; }
	return *this;
}

template<int dimension>
Point<dimension> Point<dimension>::operator*(const double s) const
{
	Point<dimension> npoint;
	for (int i = 0; i < dimension; i++){ npoint.coord[i] = coord[i] * s; }
	return npoint;
}

template<int dimension>
Point<dimension> Point<dimension>::operator*(const Point & pt) const
{
	Point<dimension> npoint;
	for (int i = 0; i < dimension; i++){ npoint.coord[i] = coord[i] * pt.coord[i]; }
	return npoint;
}

template<int dimension>
Point<dimension>& Point<dimension>::operator*=(const double s)
{
	for (int i = 0; i < dimension; i++){ coord[i] *= s; }
	return *this;
}

template<int dimension>
Point<dimension>& Point<dimension>::operator*=(const Point & pt)
{
	for (int i = 0; i < dimension; i++){coord[i] *= pt[i]; }
	return *this;
}

static double sign(double x)
{
	if (x>0.f) return 1.f;
	else if (x<0.f) return -1.f;
	else return 0.f;
}

template<int dimension>
Point<dimension> Point<dimension>::operator/(const double s) const
{
	Point<dimension> npoint;
	if (s == 0.f){
		printf("Zero Division \n");
		for (int i = 0; i < dimension; i++){
			npoint.coord[i] = DBL_MAX*sign(coord[i]);
		}
	}
	else{
		for (int i = 0; i < dimension; i++){ 
			npoint.coord[i] = coord[i] / s; 
		}
	}
	return npoint;
}

template<int dimension>
Point<dimension> Point<dimension>::operator/(const Point & pt) const
{
	Point<dimension> npoint;
	for (int i = 0; i < dimension; i++){
		if (pt.coord[i] == 0.f){ 
			printf("Zero Division \n");
			npoint.coord[i] = DBL_MAX*sign(coord[i]);
		}
		else{
			npoint.coord[i] = coord[i] / pt.coord[i]; 
		}

	}
	return npoint;
}

template<int dimension>
Point<dimension>& Point<dimension>::operator/=(const double s)
{
	if (s == 0.f){
		printf("Zero Division \n"); 
		for (int i = 0; i < dimension; i++){
			coord[i] = DBL_MAX*sign(coord[i]);
		}
	}
	else
	{
		for (int i = 0; i < dimension; i++){
			coord[i] /= s; 
		}
	}
	return *this;
}

template<int dimension>
Point<dimension>& Point<dimension>::operator/=(const Point & pt)
{
	for (int i = 0; i < dimension; i++){
		if (pt.coord[i] == 0.f){
			printf("Zero Division \n"); 
			coord[i] = DBL_MAX*sign(coord[i]);
		}
		else
		{
			coord[i] /= pt[i]; }
		}
		
	return *this;
}


template<int dimension>
Point<dimension>& Point<dimension>::normalize()
{
	if (norm()>0.f)
	{	
		(*this) *= (1.f / norm());
	}
	return *this;
}
template<int dimension>
Point<dimension> Point<dimension>::normalvector() const
{
	Point<dimension> npoint;
	if (norm() > 0.f)
	{
		npoint = (*this)*(1.f / norm());
	}
	return npoint;
}
template<int dimension>
double Point<dimension>::norm(void) const
{
	double d_norm = 0.f;
	for (int i = 0; i < dimension; i++){ d_norm += (coord[i] * coord[i]); }
	d_norm = sqrt(d_norm);
	return d_norm;
}

template<int dimension>
double Point<dimension>::linf_norm(void) const
{
	double linf_norm_val = 0;
	for (int i = 0; i < dimension; i++){
		if (abs(coord[i]) > linf_norm_val){
			linf_norm_val = abs(coord[i]);
		}
	}
	return linf_norm_val;
}

template<int dimension>
double Point<dimension>::l1_norm(void) const
{
	double l1_norm_val = 0;
	for (int i = 0; i < dimension; i++){
			l1_norm_val += abs(coord[i]);
	}
	return l1_norm_val;
}

template<int dimension>
double Point<dimension>::squared_norm(void) const
{
	double d_norm = 0.f;
	for (int i = 0; i < dimension; i++){ d_norm += (coord[i] * coord[i]); }
	return d_norm;
}

template<int dimension>
double Point<dimension>::dot(const Point<dimension> pt) const
{
	double d_dot = 0.f;
	for (int i = 0; i < dimension; i++){ d_dot += (coord[i] * pt.coord[i]); }
	return d_dot;
}

template<int dimension>
Point<dimension> Point<dimension>::cross(const Point<dimension> & pt) const
{
	Point<dimension> npoint;
	if (dimension == 3)
	{
		npoint.coord[0] = coord[1] * pt.coord[2] - coord[2] * pt.coord[1];
		npoint.coord[1] = coord[2] * pt.coord[0] - coord[0] * pt.coord[2];
		npoint.coord[2] = coord[0] * pt.coord[1] - coord[1] * pt.coord[0];
	}
	if (dimension == 2)
	{
		npoint.coord[0] = abs(coord[0] * pt.coord[1] - coord[1] * pt.coord[0]);// Area term
		npoint.coord[1] = 0;
	}
	return npoint;
}

template<int dimension>
double Point<dimension>::min_value( int & index ) const
{
	double temp_min_value = DBL_MAX;
	int temp_index = -1;
	for (int i = 0; i < dimension; i++)
	{
		if (temp_min_value > coord[i])
		{
			temp_min_value = coord[i];
			temp_index = i;
		}
	}
	index = temp_index;
	return temp_min_value;
}

template<int dimension>
double Point<dimension>::min_positive_value(int & index) const
{
	double temp_min_value = DBL_MAX;
	int temp_index = -1;
	for (int i = 0; i < dimension; i++)
	{
		if (temp_min_value > coord[i] && coord[i] > 0.f)
		{
			temp_min_value = coord[i];
			temp_index = i;
		}
	}
	index = temp_index;
	return temp_min_value;
}

template<int dimension>
double Point<dimension>::max_value(int & index) const
{
	double temp_max_value = -DBL_MAX;
	int temp_index = -1;
	for (int i = 0; i < dimension; i++)
	{
		if (temp_max_value < coord[i])
		{
			temp_max_value = coord[i];
			temp_index = i;
		}
	}
	index = temp_index;
	return temp_max_value;
}

template<int dimension>
void Point<dimension>::squaredRoot()
{
	for (int i = 0; i < dimension; i++){
		if (coord[i] < 0.f){
			printf("Squared Root Of Negative Number!!\n");
		}
		else{
			coord[i] = sqrt(coord[i]);
		}
	}
}

template<int dimension>
Point<dimension> Point<dimension>::squaredRoot() const
{
	Point<dimension> output;
	for (int i = 0; i < dimension; i++){
		if (coord[i] < 0.f){
			printf("Squared Root Of Negative Number!!\n");
		}
		else{
			output.coord[i] = sqrt(coord[i]);
		}
	}
	return output;
}

template<int dimension>
double Point<dimension>::Mean() const
{
	if (dimension > 1){
		double cum_sum = 0.f;
		for (int i = 0; i < dimension; i++){
			cum_sum += coord[i];
		}
		cum_sum /= static_cast<double>(dimension);
		return cum_sum;
	}
	else{
		return coord[0];
	}


}

template<int dimension>
Point<dimension> Point<dimension>::projection(const Point<dimension> & pt) const
{
	if (this->squared_norm() != 0.f)
	{
		return  (*this)*(this->dot(pt) / this->squared_norm());
	}
	else
	return Point<dimension>();
}


template<int dimension>
Point<dimension> CentralizedRotation(const Point<dimension> & center, const Point<dimension> & direction, const Point<dimension> & position, const double angle)
{
	Point<dimension> centralized_position = position - center;
	Point<dimension> v1 = direction.cross(centralized_position);
	Point<dimension> v2 = direction.cross(v1);
	Point<dimension> rotated_centralized_position = centralized_position + (v1*sin(angle)) + (v2*(1.f-cos(angle)));

	Point<dimension> result = center + rotated_centralized_position;

	return result;
}

void FindDiagonal(Point<2> * positions, int polygon_degree, int & diagonal_base_vertex, int & diagonal_cyclce_lenght);
bool InteriorVertex(int polygon_degree, Point<2> * positions);
int FindIndex(int input_index, int list_size, int * list_values);
bool SegmentPairIntersection(Point<2> central_position, Point<2> pivot, Point<2> edge_prev_vertex, Point<2> edge_next_vertex);
bool InsideBox(Point<3> position, Point<3> reference_vertex, Point<3> direction0, Point<3> direction1, Point<3> direction2, double lenght0, double lenght1, double lenght2);

template<int dimension>
void vectorPlanarRotation(Point<dimension>in_x, Point<dimension> in_y, Point<dimension>& out_x, Point<dimension> & out_y, Point<2> initial_direction, Point<2> final_direction)
{
	Point<2> initial_orthogonal = Point<2>(-initial_direction[1], initial_direction[0]);
	Point<2> final_orthogonal = Point<2>(-final_direction[1], final_direction[0]);

	Point<dimension> proj_principal = in_x*initial_direction[0] + in_y*initial_direction[1];
	Point<dimension> proj_orthogonal = in_x*initial_orthogonal[0] + in_y*initial_orthogonal[1];

	out_x = proj_principal * final_direction[0] + proj_orthogonal *final_orthogonal[0];
	out_y = proj_principal * final_direction[1] + proj_orthogonal *final_orthogonal[1];
}
//enum VectorMode{ COLUMN, ROW };


template< int cols,int rows>
class Matrix{
public:
	double m[cols][rows];
	Matrix(void);
	//Matrix(const Point<3>& p1, const Point<3>& p2, const Point<3>& p3, VectorMode mode);
	
	void InsertColumn(const Point<rows> p, int index);
	void InsertRow(const Point<cols> p, int index);
	static Matrix ColumnMatrix(const Point<rows> p1, ...);
	static Matrix RowMatrix(const Point<cols> p1, ...);
	Point<cols> ExtractColumn(const int index);
	Point<rows> ExtractRow(const int index);
	Point<rows> ReduceColumns();
	Point<cols> ReduceRows();

	/** This method returns the entry of the matrix in the col-th column and the row-th row.*/
	double& operator() (int col, int row);
	/** This method returns the entry of the matrix in the col-th column and the row-th row.*/
	double& index(int col, int row);

	/** This method returns the determinant of the matrix.*/
	double det(void) const;

	template<int ncols, int  nrows> Matrix<ncols, rows> LeftMultiply(const Matrix<ncols, nrows> n);
	template<int ncols, int  nrows> Matrix<cols,  nrows>  RightMultiply(const Matrix<ncols, nrows> n);
	//Point<3> Extract(VectorMode mode,const int index);

	/** This method returns the negative of the matrix */
	Matrix operator- (void) const;

	/** This method multiplies two matrices and returns the product.*/
	Matrix  operator* (const Matrix& m) const;
	/** This method multiplies the current matrix (on the right) by the input matrix..*/
	Matrix& operator*=(const Matrix& m);

	/** This method adds two matrices and returns the sum. */
	Matrix  operator+ (const Matrix& m) const;
	/** This method adds the input matrix to the current matrix. */
	Matrix& operator+=(const Matrix& m);

	/** This method subtracts two matrices and returns the difference. */
	Matrix  operator- (const Matrix& m) const;
	/** This method subtracts the input matrix from the current matrix. */
	Matrix& operator-=(const Matrix& m);

	/** This method scales the entries of a matrix and returns a new matrix. */
	Matrix  operator* (double f) const;
	/** This method scales the entries of the current matrix. */
	Matrix& operator*=(double f);
	/** This method divides the entries of a matrix and returns a new matrix. */
	Matrix  operator/ (double f) const;
	/** This method divides the entries of the current matrix. */
	Matrix& operator/=(double f);

	/** This method returns the transpose of a matrix. (Note that it does not change the entries of the current matrix.)*/
	Matrix<rows,cols> transpose(void) const;
	void transpose(const Matrix<rows, cols> n);
	//static Matrix<rows, cols> Transpose(const Matrix in);
	//template<int ncols, int nrows> Matrix<ncols, nrows> generalTranspose() const;

	/** This method returns the inverse of a matrix. (Note that it does not change the entries of the current matrix.)*/
	Matrix invert(void) const;

	/** This static method tries to invert the input matrix and write it out into the output. A value of 0 is returned if the matrix has 0 determinant.*/
	static int Invert(const Matrix& in, Matrix& out);

	/** This method transforms a 3D point. */
	Point<rows> operator*(const Point<cols>& p) const;

	/** This static method returns the identity matrix. */
	static Matrix IdentityMatrix(void);

	/** This method returns the sum of the squares of the matrix entries */
	double squareNorm(void) const;
	double frobeniusNorm(void) const;
	/** This method returns sum of the squares of the entries of the difference matrix. */
	static double SquareL2Difference(const Matrix& m1, const Matrix& m2);

	/** This method returns the exponent of a matrix */
	static Matrix Exp(const Matrix& m, int iter = 100);
	/** This method returns the logarithm of a matrix */
	static Matrix Log(const Matrix& m, double eps = 0.0001);
	/** This method returns the square-root of a matrix */
	static Matrix SquareRoot(const Matrix& m, double eps = 0.000001);

	/** This method computes the SVD decomposition of the upper 3x3 matrix, such that
	* r1 and r2 are rotations and the upper 3x3 matrix is equal to r1*diagonal*r2 */
	//void SVD(Matrix& r1, Matrix& diagonal, Matrix& r2) const;

	///** This method factors a matrix as the product of a rotation and a symmetric matrix */
	//void Factor(Matrix& rot, Matrix& sym) const;

	///** This method returns the closest 3x3 rotation matrix */
	//Matrix closestRotation(void) const;

	///** This method returns nearest symmetric matrix */
	//Matrix symmetrize(void) const;
	///** This method returns nearest skew-symmetric matrix */
	//Matrix skewSymmetrize(void) const;
};

typedef Matrix<3, 3> Matrix3D;
typedef Matrix<2, 2> Matrix2D;

//Matrix3D RankOneMatrix(Point<3> u, Point<3> v);
template<int dimension>
Matrix<3, dimension> RankOneMatrix(Point<dimension> u, Point<3> v);

Matrix3D ApplyNeighbourTransformation(Matrix3D transformation_to_neighbour, Matrix3D input, Point<3> normal, bool apply_transformation);
Point<3> ApplyNeighbourTransformation(Matrix3D transformation_to_neighbour, Point<3> input, Point<3> normal, bool apply_transformation);

void Matrix2DPrincipalComponents(Matrix2D matrix, Point<2> * principal_directions, double * eigenvalues);
void Matrix2DPrincipalComponents(Matrix2D matrix, Point<2> & principal_direction, double & principal_eigenvalue, double & secondary_eigenvalue);
#include "matrix.inl"

#endif // GEOMETRY_INCLUDED