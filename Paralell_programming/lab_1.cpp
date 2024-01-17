#include <iostream>
#include <memory>
#include <omp.h>


double average(const double* vector, size_t n)
{
	double sum = 0;
	for (size_t index = 0; index < n; index++)
	{
		sum += vector[index];
	}
	return sum / n;
}


double average_reduce(const double* vector, size_t n)
{
	double sum = 0.0;
#pragma omp parallel for reduction(+ :sum)
	for (size_t index = 0; index < n; index++)
	{
		sum += vector[index];
	}
	return sum / n;
}