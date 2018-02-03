#include "file-io.h"
#include <cstdio>
#include <stdlib.h>

void FileIO::readInputdd(char * input, double *A, int num_entries)
{

	printf("Reading From: %s \n", input);
	FILE * file;
	file = fopen(input, "r");
	if(file == 0){
		printf("File Not Found!!");
	}
	double value;
	for (int i = 0; i < num_entries; i++)
	{
			fscanf(file, "%lf", &value);
			A[i] = value;
    }
	fclose(file);
}

void FileIO::readInputff(char * input, float *A, int num_entries)
{

	printf("Reading From: %s \n", input);
	FILE * file;
	file = fopen(input, "r");
	if(file == 0){
		printf("File Not Found!!");
	}
	float value;
	for (int i = 0; i < num_entries; i++)
	{
			fscanf(file, "%f", &value);
			A[i] = value;
    }
	fclose(file);
}

void FileIO::readInputdf(char * input, float *A, int num_entries)
{

	printf("Reading From: %s \n", input);
	FILE * file;
	file = fopen(input, "r");
	if(file == 0){
		printf("File Not Found!!");
	}
    double value;
	for (int i = 0; i < num_entries; i++)
	{
			fscanf(file, "%lf", &value);
			A[i] = static_cast<float>(value);
    }
	fclose(file);
}

void FileIO::readInputii(char * input, int *A, int num_entries)
{

	printf("Reading From: %s \n", input);
	FILE * file;
	file = fopen(input, "r");
	if(file == 0){
		printf("File Not Found!!");
	}
	int value;
	for (int i = 0; i < num_entries; i++)
	{
			fscanf(file, "%d", &value);
			A[i] = value;
    }
	fclose(file);
}

void FileIO::writeOutputdd(char * output, double *A, int num_entries)
{

	printf("Writing To: %s \n", output);
	FILE * file;
	file = fopen(output, "w");
	for (int i = 0; i < num_entries; i++)
	{
			fprintf(file, "%lf ", A[i]);
	}
	fclose(file);
}

void FileIO::writeOutputff(char * output, float *A, int num_entries)
{

	printf("Writing To: %s \n", output);
	FILE * file;
	file = fopen(output, "w");
	for (int i = 0; i < num_entries; i++)
	{
			fprintf(file, "%f ", A[i]);
	}
	fclose(file);
}

void FileIO::writeOutputii(char * output, int *A, int num_entries)
{

	printf("Writing To: %s \n", output);
	FILE * file;
	file = fopen(output, "w");
	for (int i = 0; i < num_entries; i++)
	{
			fprintf(file, "%d", A[i]);
	}
	fclose(file);
} 

void FileIO::writeOutputui(char * output, unsigned char *A, int num_entries)
{

	printf("Writing To: %s \n", output);
	FILE * file;
	file = fopen(output, "w");
	for (int i = 0; i < num_entries; i++)
	{
			fprintf(file, "%d ", A[i]);
	}
	fclose(file);
}