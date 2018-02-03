#ifndef ARRAY_STATISTICS_INCLUDED
#define ARRAY_STATISTICS_INCLUDED

#include "geometry.h"

template<int dimension>
class Statistics_Array
{
public:
	Statistics_Array(Point<dimension> * p_entries, int p_array_size, int p_num_percentiles = 10){
		num_percentiles = p_num_percentiles;

		entries = p_entries;
		array_size = p_array_size;

		value_l1_norm_percentiles = new double[num_percentiles];
		value_l2_norm_percentiles = new double[num_percentiles];
		value_linf_norm_percentiles = new double[num_percentiles];

		active_percentile = num_percentiles - 2;

		verbose = false;
	}
	~Statistics_Array(){ delete[] value_l2_norm_percentiles; delete[] value_l1_norm_percentiles; delete[] value_linf_norm_percentiles; }
	int num_percentiles;
	int active_percentile;
	double max_value_l2_norm;
	double mean_value_l2_norm;
	double * value_l2_norm_percentiles;

	double max_value_linf_norm;
	double mean_value_linf_norm;
	double * value_linf_norm_percentiles;

	double max_value_l1_norm;
	double mean_value_l1_norm;
	double * value_l1_norm_percentiles;

	Point<dimension> mean_value;

	Point<dimension> * entries;
	int array_size;

	void updateStatistics();
	void setArrayStatisticsMean();
	void setArrayStatisticsL1();
	void setArrayStatisticsL2();
	void setArrayStatisticsLInf();

	bool verbose;
};


template<int dimension>
void Statistics_Array<dimension>::setArrayStatisticsL1()
{
	double  * norms = new double[array_size];
	for (int i = 0; i < array_size; i++){
		norms[i] = entries[i].l1_norm();
	}

	int sample_size = 10000;
	double  * sampled_norms = new double[sample_size];
	Sort::SampleSet(norms, array_size, sampled_norms, sample_size);
	qsort(sampled_norms, sample_size, sizeof(double), Sort::doubleCompare);
	for (int i = 0; i < num_percentiles -1; i++)
	{
		int index = static_cast<int>(floor(static_cast<double>(sample_size - 1)*static_cast<double>(i + 1) / static_cast<double>(num_percentiles)));
		value_l1_norm_percentiles[i] = sampled_norms[index];
		if (verbose){
			printf(" L1 Norm Percentile %d of %d  = %g \n", i + 1, num_percentiles, value_l1_norm_percentiles[i]);
		}
		
	}

	max_value_l1_norm = Sort::MaxValue(norms, array_size);
	value_l1_norm_percentiles[num_percentiles - 1] = max_value_l1_norm;
	if (verbose){
		printf(" L1 Norm Percentile %d of %d  = %g \n", num_percentiles, num_percentiles, value_l1_norm_percentiles[num_percentiles - 1]);
		printf("Max L1 Norm = %g \n", max_value_l1_norm);
	}

	mean_value_l1_norm = Sort::MeanValue_Cummulative(norms, array_size);
	if (verbose){
		printf("Mean L1 Norm = %g \n", mean_value_l1_norm);
	}
	delete norms;
	delete sampled_norms;
}

template<int dimension>
void Statistics_Array<dimension>::setArrayStatisticsL2()
{
	double  * norms = new double[array_size];
	for (int i = 0; i < array_size; i++){
		norms[i] = entries[i].norm();
	}

	int sample_size = 10000;
	double  * sampled_norms = new double[sample_size];
	Sort::SampleSet(norms, array_size, sampled_norms, sample_size);
	qsort(sampled_norms, sample_size, sizeof(double), Sort::doubleCompare);
	for (int i = 0; i < num_percentiles - 1; i++)
	{
		int index = static_cast<int>(floor(static_cast<double>(sample_size - 1)*static_cast<double>(i + 1) / static_cast<double>(num_percentiles)));
		value_l2_norm_percentiles[i] = sampled_norms[index];
		if (verbose){
			printf(" L2 Norm Percentile %d of %d  = %g \n", i + 1, num_percentiles, value_l2_norm_percentiles[i]);
		}
	}

	max_value_l2_norm = Sort::MaxValue(norms, array_size);
	value_l2_norm_percentiles[num_percentiles - 1] = max_value_l2_norm;
	if (verbose){
		printf(" L2 Norm Percentile %d of %d  = %g \n", num_percentiles, num_percentiles, value_l2_norm_percentiles[num_percentiles - 1]);
		printf("Max L2 Norm = %g \n", max_value_l2_norm);
	}

	mean_value_l2_norm = Sort::MeanValue_Cummulative(norms, array_size);
	if (verbose){
		printf("Mean L2 Norm = %g \n", mean_value_l2_norm);
	}
	delete norms;
	delete sampled_norms;
}

template<int dimension>
void Statistics_Array<dimension>::setArrayStatisticsLInf()
{
	double  * norms = new double[array_size];
	for (int i = 0; i < array_size; i++){
		norms[i] = entries[i].linf_norm();
	}

	int sample_size = 10000;
	double  * sampled_norms = new double[sample_size];
	Sort::SampleSet(norms, array_size, sampled_norms, sample_size);
	qsort(sampled_norms, sample_size, sizeof(double), Sort::doubleCompare);
	for (int i = 0; i < num_percentiles-1; i++)
	{
		int index = static_cast<int>(floor(static_cast<double>(sample_size - 1)*static_cast<double>(i + 1) / static_cast<double>(num_percentiles)));
		value_linf_norm_percentiles[i] = sampled_norms[index];
		if (verbose){
			printf(" LInf Norm Percentile %d of %d  = %g \n", i + 1, num_percentiles, value_linf_norm_percentiles[i]);
		}
	}

	max_value_linf_norm = Sort::MaxValue(norms, array_size);
	value_linf_norm_percentiles[num_percentiles - 1] = max_value_linf_norm;
	if (verbose){
		printf(" LInf Norm Percentile %d of %d  = %g \n", num_percentiles, num_percentiles, value_linf_norm_percentiles[num_percentiles - 1]);
		printf("Max LInf Norm = %g \n", max_value_linf_norm);
	}

	mean_value_linf_norm = Sort::MeanValue_Cummulative(norms, array_size);
	if (verbose){
		printf("Mean LInf Norm = %g \n", mean_value_linf_norm);
	}
	delete norms;
	delete sampled_norms;
}

template<int dimension>
void Statistics_Array<dimension>::setArrayStatisticsMean()
{
	Point<dimension> cummulative_mean;
	for (int i = 0; i < array_size; i++){
		cummulative_mean += entries[i];
	}
	mean_value = cummulative_mean / static_cast<double>(array_size);
}


template<int dimension>
void Statistics_Array<dimension>::updateStatistics()
{
	setArrayStatisticsMean();
	setArrayStatisticsL1();
	setArrayStatisticsL2();
	setArrayStatisticsLInf();
}

#endif // ARRAY_STATISTICS_INCLUDED