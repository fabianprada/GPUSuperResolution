/*
Copyright (c) 2018, Fabian Prada
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution.

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/


#ifndef CPU_SOLVER_INCLUDE
#define CPU_SOLVER_INCLUDE

#include <vector>
#include "Eigen\Sparse"


typedef Eigen::Triplet<double> Triplet_d;
typedef std::vector<Triplet_d> Triplet_dList;
typedef Eigen::SparseMatrix<double> SparseMatrix_d; // declares a column-major sparse matrix type of double
typedef Eigen::SimplicialCholesky<SparseMatrix_d> CholeskyFactorization_d;
typedef Eigen::VectorXd Vector_d;

typedef Eigen::Triplet<float> Triplet_f;
typedef std::vector<Triplet_f> Triplet_fList;
typedef Eigen::SparseMatrix<float> SparseMatrix_f; // declares a column-major sparse matrix type of float
typedef Eigen::SimplicialCholesky<SparseMatrix_f> CholeskyFactorization_f;
typedef Eigen::VectorXf Vector_f;

namespace EigenSolver{

	CholeskyFactorization_d * SetLaplaceCholeskyFactorization_d(const unsigned int array_width,const unsigned int array_height);
	CholeskyFactorization_d * SetCholeskyFactorization_d(double  * filter, const int filter_width,const int filter_height, const unsigned int array_width,const unsigned int array_height);
	SparseMatrix_d  GenerateSparseMatrix_d(double  * filter, const int filter_width,const int filter_height, const unsigned int array_width,const unsigned int array_height);
	void SolveFromCholesky_d(CholeskyFactorization_d * factorization, double * buffer,const unsigned int buffer_size,  double  input_dc);
	Vector_d GenerateVector_d(double * buffer, const unsigned int buffer_size);
	void CopyFromVector_d(double * buffer, const Vector_d & values , const unsigned int buffer_size,  double scalar_correction);

	CholeskyFactorization_f * SetLaplaceCholeskyFactorization_f(const unsigned int array_width,const unsigned int array_height);
	CholeskyFactorization_f * SetCholeskyFactorization_f(float  * filter, const int filter_width,const int filter_height, const unsigned int array_width,const unsigned int array_height);
	SparseMatrix_f  GenerateSparseMatrix_f(float  * filter, const int filter_width,const int filter_height, const unsigned int array_width,const unsigned int array_height);
	void SolveFromCholesky_f(CholeskyFactorization_f * factorization, float * buffer,const unsigned int buffer_size,  float  input_dc);
	Vector_f GenerateVector_f(float * buffer, const unsigned int buffer_size);
	void CopyFromVector_f(float * buffer, const Vector_f & values , const unsigned int buffer_size,  float scalar_correction);

}

#endif // CPU_SOLVER_INCLUDE