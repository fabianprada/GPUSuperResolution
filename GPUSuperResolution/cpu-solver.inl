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
#include <stdlib.h>

template<int dimension>
void SparseSystem<dimension>::AssignNeigbourRelation(int index1, int index2, double correlation_coefficient)
{
	elements[index1].neighbours.push_back(NeighbourRelation(correlation_coefficient,index2));
	elements[index2].neighbours.push_back(NeighbourRelation(correlation_coefficient,index1));
}

template<int dimension>
void SparseSystem<dimension>::EvaluateSquaredResidue(double * x)
{
	double  * Ax = new double[num_elements];
	double  * b = new double[num_elements];
	double  * r = new double[num_elements];

	for (int channel = 0; channel < dimension; channel++)
	{
		int index = 0;
		for (const_element_iterator element_iter = elements.begin(); element_iter != elements.end(); element_iter++){
			const SystemElement<dimension> & element = *element_iter;
			b[index] = element.constraint_coefficient[channel];
			index++;
		}
		ApplyLinearOperator(x, Ax);

		double squared_residue = 0.f;
		//#pragma omp parallel for num_threads( threads ) reduction ( + : delta_new )
		for (int i = 0; i < num_elements; i++){
			r[i] = b[i] - Ax[i]; 
			squared_residue += r[i] * r[i];
		}
		printf("Squared Residue Channel %d = %f \n",channel, squared_residue);
	}
	delete[] b, delete[] d; delete[] r; delete[] Ax, delete[] Ad; delete[] x;
}

template<int dimension>
void SparseSystem<dimension>::ApplyLinearOperator(double * input, double * output)
{
	int index = 0;
	for (const_element_iterator element_iter = elements.begin(); element_iter != elements.end(); element_iter++)
	{
		const SystemElement<dimension> & element = *element_iter;
		double cum_result = input[index] * element.central_coefficient;
		for (const_neighbour_iterator neighbour_iter = element.neighbours.begin(); neighbour_iter != element.neighbours.end(); neighbour_iter++)
		{
			const NeighbourRelation & neighbour = *neighbour_iter;
			cum_result += (neighbour.neighbour_coefficient*input[neighbour.neighbour_index]);
		}
		output[index] = cum_result;
		index++;
	}
}

template<int dimension>
double  * SparseSystem<dimension>::CGSolver(double eps, int iters, bool correct_dc, bool clamp)
{
	double  * solution = new double[num_elements*dimension];

	double  * x = new double[num_elements];
	double  * Ax = new double[num_elements];
	double  * b = new double[num_elements];
	double  * d = new double[num_elements];
	double  * Ad = new double[num_elements];
	double  * r = new double[num_elements];


	for (int channel = 0; channel < dimension; channel++)
	{

#pragma omp parallel for
		for (int i = 0; i < num_elements; i++){
			x[i] = elements[i].initial_value[channel];
			b[i] = elements[i].constraint_coefficient[channel];
		}

		double initial_dc = 0.f;
		if (correct_dc)
		{
#pragma omp parallel for reduction ( + : initial_dc)
			for (int i = 0; i < num_elements; i++){ 
				initial_dc += x[i]; 
			} //Implement Binary sum
			printf("Initial DC Channel %d = %g \n", channel, initial_dc);
		}


		double delta_new = 0.f;

		ApplyLinearOperator(x, Ax);

#pragma omp parallel for reduction ( + : delta_new )
		for (int i = 0; i < num_elements; i++){
			d[i] = r[i] = b[i] - Ax[i];
			delta_new += r[i] * r[i];
		}
		printf("Initial Residue Channel %d = %g \n", channel, delta_new);
		//if (delta_new<eps)
		//{
		//	delete[] r, delete[] d, delete[] b; delete[] x;
		//	return x;
		//}
		int ii = 0;
		if (delta_new > eps)
		{
			double delta_0 = delta_new;
			for (ii = 0; ii<iters && delta_new>eps*delta_0; ii++)
			{
				//printf("Squared Residue Channel %d Iteration %d  = %f \n", channel, ii, delta_new);
				//L(d, q);
				ApplyLinearOperator(d, Ad);
				double dDotAd = 0;
#pragma omp parallel for reduction( + : dDotAd )
				for (int i = 0; i < num_elements; i++){ 
					dDotAd += d[i] * Ad[i]; 
				}
				double alpha = delta_new / dDotAd;
#pragma omp parallel for
				for (int i = 0; i < num_elements; i++){ 
					x[i] += d[i] * alpha; 
				}
				double delta_old = delta_new;
				delta_new = 0;

				const int RESET_COUNT = 50;
				if ((ii%RESET_COUNT) == (RESET_COUNT - 1))
				{
					ApplyLinearOperator(x, Ax);
					for (int i = 0; i < num_elements; i++){ d[i] = r[i] = b[i] - Ax[i]; delta_new += r[i] * r[i]; }
				}
				else
				{
#pragma omp parallel for reduction( + : delta_new )
					for (int i = 0; i < num_elements; i++){
						r[i] -= Ad[i] * alpha;   delta_new += r[i] * r[i]; 
					}
					double beta = delta_new / delta_old;
#pragma omp parallel for 
					for (int i = 0; i < num_elements; i++){
						d[i] = r[i] + d[i] * beta; 
					}
				}
			}
		}
		printf("Num Iterations Channel %d = %d \n", channel, ii);
		printf("Final Residue Channel %d = %g \n", channel, delta_new);

		if (clamp)
		{
			printf("Clamping Solution \n");
#pragma omp parallel for 
			for (int i = 0; i < num_elements; i++){ 
				x[i] = x[i] > 1.f ? 1.f : x[i]; x[i] = x[i] < 0.f ? 0.f : x[i]; 
			}
		}

		double final_dc = 0.f;
		if (correct_dc)
		{
#pragma omp parallel for reduction ( + : final_dc )
			for (int i = 0; i < num_elements; i++){ 
				final_dc += x[i]; 
			}
			printf("Final DC Channel %d = %g \n", channel, final_dc);
			printf("Correcting DC \n");
			double dc_correction = (initial_dc - final_dc) / static_cast<double>(num_elements);
#pragma omp parallel for
			for (int i = 0; i < num_elements; i++){ 
				x[i] += dc_correction; 
			}
		}
#pragma omp parallel for
		for (int i = 0; i < num_elements; i++){
			solution[i*dimension + channel] = x[i];
		}
		
	}
	delete[] b, delete[] d; delete[] r; delete[] Ax, delete[] Ad; delete[] x;

	return solution;
}

//template<int dimension>
//double  * SparseSystem<dimension>::CGSolver(double eps, int iters, bool correct_dc, bool clamp)
//{
//	double  * solution = new double[num_elements*dimension];
//
//	double  * x = new double[num_elements];
//	double  * Ax = new double[num_elements];
//	double  * b = new double[num_elements];
//	double  * d = new double[num_elements];
//	double  * Ad = new double[num_elements];
//	double  * r = new double[num_elements];
//
//
//	for (int channel = 0; channel < dimension; channel++)
//	{
//		int index = 0;
//		for (const_element_iterator element_iter = elements.begin(); element_iter != elements.end(); element_iter++){
//			const SystemElement<dimension> & element = *element_iter;
//			x[index] = element.initial_value[channel];
//			b[index] = element.constraint_coefficient[channel];
//			index++;
//		}
//
//		double initial_dc = 0.f;
//		if (correct_dc)
//		{
//			for (int i = 0; i < num_elements; i++){ initial_dc += x[i]; } //Implement Binary sum
//			printf("Initial DC Channel %d = %f \n", channel, initial_dc);
//		}
//
//
//		double delta_new = 0.f;
//
//		ApplyLinearOperator(x, Ax);
//
//		//#pragma omp parallel for num_threads( threads ) reduction ( + : delta_new )
//		for (int i = 0; i < num_elements; i++){ d[i] = r[i] = b[i] - Ax[i]; delta_new += r[i] * r[i]; }
//
//		//if (delta_new<eps)
//		//{
//		//	delete[] r, delete[] d, delete[] b; delete[] x;
//		//	return x;
//		//}
//
//		double delta_0 = delta_new;
//		int ii;
//		for (ii = 0; ii<iters && delta_new>eps*delta_0; ii++)
//		{
//			printf("Squared Residue Channel %d = %f \n", channel, delta_new);
//			//L(d, q);
//			ApplyLinearOperator(d, Ad);
//			double dDotAd = 0;
//			//#pragma omp parallel for num_threads( threads ) reduction( + : dDotQ )
//			for (int i = 0; i < num_elements; i++){ dDotAd += d[i] * Ad[i]; }
//			double alpha = delta_new / dDotAd;
//			//#pragma omp parallel for num_threads( threads )
//			for (int i = 0; i < num_elements; i++){ x[i] += d[i] * alpha; }
//			double delta_old = delta_new;
//			delta_new = 0;
//
//			const int RESET_COUNT = 50;
//			if ((ii%RESET_COUNT) == (RESET_COUNT - 1))
//			{
//				ApplyLinearOperator(x, Ax);;
//				//#pragma omp parallel for num_threads( threads ) reduction ( + : delta_new )
//				for (int i = 0; i < num_elements; i++){ d[i] = r[i] = b[i] - Ax[i];   delta_new += r[i] * r[i]; }
//			}
//			else
//			{
//				//#pragma omp parallel for num_threads( threads ) reduction( + : delta_new )
//				for (int i = 0; i < num_elements; i++){ r[i] -= Ad[i] * alpha;   delta_new += r[i] * r[i];  x[i] += d[i] * alpha; }
//			}
//			double beta = delta_new / delta_old;
//			//#pragma omp parallel for num_threads( threads )
//			for (int i = 0; i<num_elements; i++){ d[i] = r[i] + d[i] * beta; }
//
//		}
//		printf("Num Iterations Channel %d = %d \n", channel,ii);
//		if (clamp)
//		{
//			for (int i = 0; i < num_elements; i++){ x[i] = x[i] > 1.f ? 1.f : x[i]; x[i] = x[i] < 0.f ? 0.f : x[i]; }
//		}
//
//		double final_dc = 0.f;
//		if (correct_dc)
//		{
//			for (int i = 0; i < num_elements; i++){ final_dc += x[i]; }
//			printf("Final DC Channel %d = %f \n", channel, final_dc);
//			double dc_correction = (initial_dc - final_dc) / static_cast<double>(num_elements);
//			for (int i = 0; i < num_elements; i++){ x[i] += dc_correction; }
//		}
//		for (int i = 0; i < num_elements; i++){ 
//			solution[i*dimension+channel] = x[i]; 
//		}
//		printf("Final Residue Channel %d = %f \n", channel, delta_new);
//	}
//	delete[] b, delete[] d; delete[] r; delete[] Ax, delete[] Ad; delete[] x;
//
//	return solution;
//}


template<int dimension>
SparseMatrixd SparseSystem<dimension>::GenerateSparseMatrix()
{
	TripletdList triplet_list;
	int num_coefficients = elements.size();
	for (int i = 0; i < elements.size(); i++){
		num_coefficients += elements[i].neighbours.size();
	}

	triplet_list.reserve(num_coefficients);
	for (int i = 0; i < elements.size(); i++){
		SystemElement<dimension> & current_element = elements[i];
		Neighbourhood & neighbours = current_element.neighbours;
		triplet_list.push_back(Tripletd(current_element.element_index, current_element.element_index, current_element.central_coefficient));
		for (int j = 0; j < neighbours.size(); j++){
			triplet_list.push_back(Tripletd(current_element.element_index, neighbours[j].neighbour_index, neighbours[j].neighbour_coefficient));
		}
	}

	SparseMatrixd A(elements.size(), elements.size());
	A.setFromTriplets(triplet_list.begin(), triplet_list.end());
	return A;
}

template<int dimension>
Vectord SparseSystem<dimension>::GenerateRHS(int channel)
{
	Vectord b(elements.size());
	for (int i = 0; i < elements.size(); i++){
		b(i) = elements[i].constraint_coefficient[channel];
	}
	return b;
}

template<int dimension>
double SparseSystem<dimension>::InitialDC(int channel)
{
	double dc = 0.f;
	for (int i = 0; i < elements.size(); i++){
		dc += elements[i].initial_value[channel];
	}
	return dc;
}