#include "geometry.h"

template< int cols, int rows>
Matrix<cols, rows>::Matrix(void){
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			m[i][j] = 0.f;
		}
	}
}

template< int cols, int rows>
void Matrix<cols, rows>::InsertColumn(const Point<rows> p, int index)
{
	if (index < 0 || index>cols - 1){ printf("Column Index Out of Bounds \n"); }
	for (int j = 0; j < rows; j++){
		m[index][j] = p[j];
	}
}

template<int dimension>
Matrix<3, dimension> RankOneMatrix(Point<dimension> u, Point<3> v)
{
	Matrix<3, dimension>  matrix;
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			matrix.m[j][i] = u[i] * v[j];
		}
	}
	return matrix;
}

template< int cols, int rows>
void Matrix<cols, rows>::InsertRow(const Point<cols> p, int index)
{
	if (index < 0 || index>rows - 1){ printf("Row Index Out of Bounds \n"); }
	for (int i = 0; i < cols; i++){
		m[i][index] = p[i];
	}
}

template< int cols, int rows>
Matrix<cols, rows> Matrix<cols, rows>::ColumnMatrix(const Point<rows> p1, ...)
{
	Matrix<cols,rows> o;
	o.InsertColumn(p1, 0);
	va_list vl;
	va_start(vl, p1);
	for (int i = 1; i<cols; i++) o.InsertColumn(va_arg(vl, Point<rows>), i);
	va_end(vl);
	return o;
}

template< int cols, int rows>
Matrix<cols, rows> Matrix<cols, rows>::RowMatrix(const Point<cols> p1, ...)
{
	Matrix<cols, rows> o;
	o.InsertRow(p1, 0);
	va_list vl;
	va_start(vl, p1);
	for (int i = 1; i<rows; i++) o.InsertRow(va_arg(vl, Point<cols>), i);
	va_end(vl);
	return o;
}


template< int cols, int rows>
Point<cols> Matrix<cols, rows>::ExtractColumn(const int index)
{
	if (index<0 || index>cols - 1){ printf("Out of bound column index \n"); }
	Point<rows> column;
	for (int j = 0; j < rows; j++) column.coord[j] = m[index][j];
	return column;
}


template< int cols, int rows>
Point<rows> Matrix<cols, rows>::ReduceColumns()
{
	Point<rows> column;
	for (int j = 0; j < rows; j++){
		double val = 0.f;
		for (int i = 0; i < cols; i++){
			val += m[i][j];
		}
		column.coord[j] = val;
	}
	return column;
}

template< int cols, int rows>
Point<cols> Matrix<cols, rows>::ReduceRows()
{
	Point<cols> row;
		for (int i = 0; i < cols; i++){
		double val = 0.f;
			for (int j = 0; j < rows; j++){
			val += m[i][j];
		}
		row.coord[i] = val;
	}
	return row;
}

template< int cols, int rows>
Point<rows> Matrix<cols, rows>::ExtractRow(const int index)
{
	if (index<0 || index>rows - 1){ printf("Out of bound column index \n"); }
	Point<cols> row;
	for (int i = 0; i < cols; i++) row.coord[i] = m[i][index];
	return row;
}

template< int cols, int rows>
double& Matrix<cols, rows>::index(int i, int j){ return m[i][j]; }

template< int cols, int rows>
double& Matrix<cols, rows>::operator() (int i, int j){ return index(i, j); }

template< int cols, int rows>
double  Matrix<cols, rows>::det(void) const {
	if (cols == 3 && rows == 3){
		return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[1][0] * (m[0][1] * m[2][2] - m[0][2] * m[2][1]) + m[2][0] * (m[0][1] * m[1][2] - m[0][2] * m[1][1]);
	}
	else
	{
		printf("Non defined determinant \n");
		return 0;
	}
}

template< int cols, int rows>
Matrix<cols, rows> Matrix<cols, rows>::operator - (void) const {
	Matrix<cols, rows> n;
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			n.m[i][j] = -m[i][j];
		}
	}
	return n;
}

template< int cols, int rows>
Matrix<cols, rows>  Matrix<cols, rows>::operator* (const Matrix<cols, rows>& n) const {
	Matrix<cols, rows> o;
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			double val = 0.f;
			for (int k = 0; k < rows; k++){
				val += m[k][j] * n.m[i][k];
			}
			o.m[i][j] = val;
		}
	}
	return o;
}

template< int cols, int rows>
template<int ncols, int nrows>
Matrix<ncols, rows> Matrix<cols, rows>::LeftMultiply(const Matrix<ncols, nrows> n)
{
	Matrix<ncols, rows> o;
	if (cols != nrows)
	{
		printf("Matrices dimensons do not match \n");
	}
	for (int i = 0; i < ncols; i++){
		for (int j = 0; j < rows; j++){
			double val = 0.f;
			for (int k = 0; k < cols; k++){
				val += m[k][j] * n.m[i][k];
			}
			o.m[i][j] = val;
		}
	}
	return o;
}


template< int cols, int rows>
template<int ncols, int nrows>
Matrix<cols, nrows> Matrix<cols, rows>::RightMultiply(const Matrix<ncols, nrows> n)
{
	Matrix<cols, nrows> o;
	if (ncols != rows)
	{
		printf("Matrices dimensons do not match \n");
	}
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < nrows; j++){
			double val = 0.f;
			for (int k = 0; k < rows; k++){
				val += m[i][k] * n.m[k][j];
			}
			o.m[i][j] = val;
		}
	}
	return o;
}

//template< int cols, int rows>
//template<int ncols, int nrows>
//Matrix<ncols, nrows> Matrix<cols, rows>::generalTranspose() const
//{
//	if (ncols != rows || nrows != cols)
//	{
//		printf("Matrix dimensions mismatch\n");
//	}
//
//	Matrix<ncols, nrows> n = Matrix<ncols, nrows>();
//	for (int i = 0; i < cols; i++){
//		for (int j = 0; j < rows; j++){
//			n.m[j][i] = m[i][j];
//		}
//	}
//	return n;
//}

template< int cols, int rows>
Matrix<rows, cols> Matrix<cols, rows>::transpose(void) const  {
	Matrix<rows, cols> o;
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			o.m[j][i] = m[i][j];
		}
	}
	return o;
}

template< int cols, int rows>
void Matrix<cols, rows>::transpose(const Matrix<rows, cols> n){
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			m[i][j] = n.m[j][i];
		}
	}
}


template< int cols, int rows>
Matrix<cols, rows>& Matrix<cols, rows>::operator*=(const Matrix& m){
	Matrix<cols, rows> temp;
	temp = *this;
	*this = temp*m;
	return *this;
}

template< int cols, int rows>
Matrix<cols, rows>  Matrix<cols, rows>::operator+ (const Matrix& n) const {
	Matrix<rows, cols> o;
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			o.m[i][j] = m[i][j] + n.m[i][j];
		}
	}
	return o;
}

template< int cols, int rows>
Matrix<cols, rows>& Matrix<cols, rows>::operator+=(const Matrix& n){
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			m[i][j] += n.m[i][j];
		}
	}
	return *this;
}

template< int cols, int rows>
Matrix<cols, rows>  Matrix<cols, rows>::operator- (const Matrix& n) const {
	Matrix o;
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			o.m[i][j] = m[i][j] - n.m[i][j];
		}
	}
	return o;
}

template< int cols, int rows>
Matrix<cols, rows>& Matrix<cols, rows>::operator-=(const Matrix& n){
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			m[i][j] -= n.m[i][j];
		}
	}
	return *this;
}


template< int cols, int rows>
Matrix<cols, rows>  Matrix<cols, rows>::operator* (double f) const {
	Matrix o;
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			o.m[i][j] = f* m[i][j];
		}
	}
	return o;
}

template< int cols, int rows>
Matrix<cols, rows>& Matrix<cols, rows>::operator*=(double f){
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			m[i][j] *= f;
		}
	}
	return *this;
}

template< int cols, int rows>
Matrix<cols, rows>  Matrix<cols, rows>::operator/ (double f) const {
	Matrix o;
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			o.m[i][j] = m[i][j] / f;
		}
	}
	return o;
}

template< int cols, int rows>
Matrix<cols, rows>& Matrix<cols, rows>::operator/=(double f){
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			m[i][j] /= f;
		}
	}
	return *this;
}


template< int cols, int rows>
Matrix<cols, rows> Matrix<cols, rows>::invert(void) const {


	if (cols == 3 && rows == 3){

		Matrix n;
		double d;

		n.m[0][0] = m[1][1] * m[2][2] - m[1][2] * m[2][1];
		n.m[1][0] = -m[0][1] * m[2][2] + m[0][2] * m[2][1];
		n.m[2][0] = m[0][1] * m[1][2] - m[0][2] * m[1][1];
		n.m[0][1] = -m[1][0] * m[2][2] + m[1][2] * m[2][0];
		n.m[1][1] = m[0][0] * m[2][2] - m[0][2] * m[2][0];
		n.m[2][1] = -m[0][0] * m[1][2] + m[0][2] * m[1][0];
		n.m[0][2] = m[1][0] * m[2][1] - m[1][1] * m[2][0];
		n.m[1][2] = -m[0][0] * m[2][1] + m[0][1] * m[2][0];
		n.m[2][2] = m[0][0] * m[1][1] - m[0][1] * m[1][0];
		d = m[0][0] * n.m[0][0] + m[1][0] * n.m[1][0] + m[2][0] * n.m[2][0];
		if (d!=0.f){
			return (n / d).transpose();
		}
		else{
			printf("Singular Matrix \n");
			printf("%g %g %g \n", m[0][0] * m[0][0] + m[1][0] * m[1][0] + m[2][0] * m[2][0],
								 m[0][1] * m[0][1] + m[1][1] * m[1][1] + m[2][1] * m[2][1], 
								 m[0][2] * m[0][2] + m[1][2] * m[1][2] + m[2][2] * m[2][2]);
			return n;
		}
	}
	else
	{
		printf("Non defined invert \n");
		return Matrix<cols,rows>();
	}
}

template< int cols, int rows>
int Matrix<cols, rows>::Invert(const Matrix<cols, rows>& m, Matrix<cols, rows>& n){


	if (cols == 3 && rows == 3){
		double d;

		n.m[0][0] = m.m[1][1] * m.m[2][2] - m.m[1][2] * m.m[2][1];
		n.m[1][0] = -m.m[0][1] * m.m[2][2] + m.m[0][2] * m.m[2][1];
		n.m[2][0] = m.m[0][1] * m.m[1][2] - m.m[0][2] * m.m[1][1];
		n.m[0][1] = -m.m[1][0] * m.m[2][2] + m.m[1][2] * m.m[2][0];
		n.m[1][1] = m.m[0][0] * m.m[2][2] - m.m[0][2] * m.m[2][0];
		n.m[2][1] = -m.m[0][0] * m.m[1][2] + m.m[0][2] * m.m[1][0];
		n.m[0][2] = m.m[1][0] * m.m[2][1] - m.m[1][1] * m.m[2][0];
		n.m[1][2] = -m.m[0][0] * m.m[2][1] + m.m[0][1] * m.m[2][0];
		n.m[2][2] = m.m[0][0] * m.m[1][1] - m.m[0][1] * m.m[1][0];
		d = m.m[0][0] * n.m[0][0] + m.m[1][0] * n.m[1][0] + m.m[2][0] * n.m[2][0];
		if (d != 0){
			n /= d;
			n = n.transpose();
			return 1;
		}
		else{ 
			printf("Singular Matrix \n");
			return 0; 
		}
	}
	else
	{
		printf("Non defined invert \n");
		return 0;
	}
}


template< int cols, int rows>
Point<rows> Matrix<cols, rows>::operator*(const Point<cols>& p) const {
	Point<rows> product;
	for (int j = 0; j < rows; j++){
		double val = 0.f;
		for (int i = 0; i < cols; i++){
			val += p[i] * m[i][j];
		}
		product.coord[j] = val;
	}
	return product;
}


template< int cols, int rows>
Matrix<cols, rows> Matrix<cols, rows>::IdentityMatrix(void){
	if (cols != nrows)
	{
		printf("Undefined Identity");
	}
	Matrix<cols, rows> o;
	for (int i = 0; i < cols; i++){
			o.m[i][i] = 1.f;
	}
	return o;
}

template< int cols, int rows>
double Matrix<cols, rows>::squareNorm(void) const{
	double n = 0;
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			n += m[i][j] * m[i][j];
		}
	}
	return n;
}

template< int cols, int rows>
double Matrix<cols, rows>::frobeniusNorm(void) const{
	double n = this->squareNorm();
	return sqrt(n);
}

template< int cols, int rows>
double Matrix<cols, rows>::SquareL2Difference(const Matrix<cols, rows>& m1, const Matrix<cols, rows>& m2){
	double n = 0;
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			double temp = m1.m[i][j] - m2.m[i][j];
			n += temp*temp;
		}
	}
	return n;
}