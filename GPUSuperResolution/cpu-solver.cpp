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


#include "cpu-solver.h"
#include "extension.h"
#include "time.h"


SparseMatrix_d EigenSolver::GenerateSparseMatrix_d(double  * filter, const int filter_width_radius,const int filter_height_radius,const unsigned int array_width,const unsigned int array_height)
{
	Triplet_dList triplet_list;
	triplet_list.reserve(filter_width_radius*filter_height_radius*array_width*array_height);
	for (int i = 0; i < array_width; i++){
		for (int j = 0; j < array_height; j++){
		    for(int ki = -filter_width_radius ; ki<=filter_width_radius; ki++){
				for(int kj = -filter_height_radius ; kj<=filter_height_radius; kj++){

					int ni = ExtensionTools::ExtensionIndex(i+ki,array_width);
					int nj = ExtensionTools::ExtensionIndex(j+kj,array_height);
					int filter_pos = (ki+filter_width_radius) + (kj+filter_height_radius)*(2*filter_width_radius+1);
					triplet_list.push_back(Triplet_d(i + j*array_width,ni + nj*array_width,filter[filter_pos]));
				}
			}
		}
	}
	SparseMatrix_d A(array_width*array_height,array_width*array_height);
	A.setFromTriplets(triplet_list.begin(), triplet_list.end());
	return A;
}



CholeskyFactorization_d * EigenSolver::SetCholeskyFactorization_d(double  * filter, const int filter_width_radius,const int filter_height_radius,const unsigned int array_width,const unsigned int array_height)
{
	SparseMatrix_d A =  GenerateSparseMatrix_d(filter,filter_width_radius,filter_height_radius,array_width,array_height);
	double start_time = GetTime();
	CholeskyFactorization_d * cholesky_factorization = new CholeskyFactorization_d(A);
	printf("Cholesky Factorization = %.4f(s)\n", GetTime() - start_time);
	return cholesky_factorization;
}

CholeskyFactorization_d * EigenSolver::SetLaplaceCholeskyFactorization_d(const unsigned int array_width,const unsigned int array_height)
{
	const int filter_width_radius =1;
	const int filter_height_radius =1;
	double * filter_values = new double[9];
	filter_values[0] = 0.f;
	filter_values[1] = -1.f;
	filter_values[2] = 0.f;
	filter_values[3] = -1.f;
	filter_values[4] = 4.f;
	filter_values[5] = -1.f;
	filter_values[6] = 0.f;
	filter_values[7] = -1.f;
	filter_values[8] = 0.f;

	CholeskyFactorization_d * factorization = SetCholeskyFactorization_d(filter_values,filter_width_radius,filter_height_radius,array_width,array_height);

	delete filter_values;
	return factorization;
}

Vector_d EigenSolver::GenerateVector_d(double * buffer, const unsigned int buffer_size)
{
	Vector_d b(buffer_size);
	for (int i = 0; i < buffer_size; i++){
		b(i) = buffer[i];
	}
	return b;
}

void EigenSolver::CopyFromVector_d(double * buffer, const Vector_d & b , const unsigned int buffer_size, double scalar_correction)
{
	for (int i = 0; i < buffer_size; i++){
		 buffer[i] =b(i) + scalar_correction;
	}
}

void EigenSolver::SolveFromCholesky_d(CholeskyFactorization_d * factorization, double * buffer, const unsigned int buffer_size, double input_dc)
{
	Vector_d rhs_vector = GenerateVector_d(buffer,buffer_size);
	double start_time = GetTime();
	Vector_d solution = factorization->solve(rhs_vector);
	printf("Cholesky Solution = %.4f(s)\n", GetTime() - start_time);
	double final_dc = solution.sum();
	printf("Final DC %f \n", final_dc);
	double scalar_correction = (input_dc -final_dc) / static_cast<double>(buffer_size);
	CopyFromVector_d(buffer,solution,buffer_size,scalar_correction);
}


////////////////////////////////////////////////////////////////////////////////////

SparseMatrix_f EigenSolver::GenerateSparseMatrix_f(float  * filter, const int filter_width_radius,const int filter_height_radius,const unsigned int array_width,const unsigned int array_height)
{
	Triplet_fList triplet_list;
	triplet_list.reserve(filter_width_radius*filter_height_radius*array_width*array_height);
	for (int i = 0; i < array_width; i++){
		for (int j = 0; j < array_height; j++){
		    for(int ki = -filter_width_radius ; ki<=filter_width_radius; ki++){
				for(int kj = -filter_height_radius ; kj<=filter_height_radius; kj++){

					int ni = ExtensionTools::ExtensionIndex(i+ki,array_width);
					int nj = ExtensionTools::ExtensionIndex(j+kj,array_height);
					int filter_pos = (ki+filter_width_radius) + (kj+filter_height_radius)*(2*filter_width_radius+1);
					triplet_list.push_back(Triplet_f(i + j*array_width,ni + nj*array_width,filter[filter_pos]));
				}
			}
		}
	}

	SparseMatrix_f A(array_width*array_height,array_width*array_height);
	A.setFromTriplets(triplet_list.begin(), triplet_list.end());
	return A;
}



CholeskyFactorization_f * EigenSolver::SetCholeskyFactorization_f(float  * filter, const int filter_width_radius,const int filter_height_radius,const unsigned int array_width,const unsigned int array_height)
{
	SparseMatrix_f A =  GenerateSparseMatrix_f(filter,filter_width_radius,filter_height_radius,array_width,array_height);
	float start_time = GetTime();
	CholeskyFactorization_f * cholesky_factorization = new CholeskyFactorization_f(A);
	printf("Cholesky Factorization = %.4f(s)\n", GetTime() - start_time);
	return cholesky_factorization;
}

CholeskyFactorization_f * EigenSolver::SetLaplaceCholeskyFactorization_f(const unsigned int array_width,const unsigned int array_height)
{
	const int filter_width_radius =1;
	const int filter_height_radius =1;
	float * filter_values = new float[9];
	filter_values[0] = 0.f;
	filter_values[1] = -1.f;
	filter_values[2] = 0.f;
	filter_values[3] = -1.f;
	filter_values[4] = 4.f;
	filter_values[5] = -1.f;
	filter_values[6] = 0.f;
	filter_values[7] = -1.f;
	filter_values[8] = 0.f;

	CholeskyFactorization_f * factorization = SetCholeskyFactorization_f(filter_values,filter_width_radius,filter_height_radius,array_width,array_height);

	delete filter_values;
	return factorization;
}

Vector_f EigenSolver::GenerateVector_f(float * buffer, const unsigned int buffer_size)
{
	Vector_f b(buffer_size);
	for (int i = 0; i < buffer_size; i++){
		b(i) = buffer[i];
	}
	return b;
}

void EigenSolver::CopyFromVector_f(float * buffer, const Vector_f & b , const unsigned int buffer_size, float scalar_correction)
{
	for (int i = 0; i < buffer_size; i++){
		 buffer[i] =b(i) + scalar_correction;
	}
}

void EigenSolver::SolveFromCholesky_f(CholeskyFactorization_f * factorization, float * buffer, const unsigned int buffer_size, float input_dc)
{
	Vector_f rhs_vector = GenerateVector_f(buffer,buffer_size);
	printf("Initial DC %f \n", rhs_vector.sum());
	float start_time = GetTime();
	Vector_f solution = factorization->solve(rhs_vector);
	printf("Cholesky Solution = %.4f(s)\n", GetTime() - start_time);
	float final_dc = solution.sum();
	printf("Final DC %f \n", final_dc);
	float scalar_correction = (input_dc -final_dc) / static_cast<float>(buffer_size);
	CopyFromVector_f(buffer,solution,buffer_size,scalar_correction);
}