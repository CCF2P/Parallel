/*#include <iostream>
#include <memory>
#include <stdalign.h>
#include <omp.h>
#include <mutex>
#include <vector>


double average_omp(const double* vector, size_t n)
{
    double sum = 0.0;
    #pragma omp parallel
    {
        unsigned t = omp_get_thread_num();
        unsigned T = omp_get_num_threads();
        for (size_t index = t; index < n; index += T)
        {
            sum += vector[index];
        }
    }
    return sum / n;
}

double average_omp_with_array(const double* vector, size_t n)
{
    double sum = 0.0;
    double* partical_sums = (double*)calloc(omp_get_num_procs(), sizeof(double));
    #pragma omp parallel
    {
        unsigned t = omp_get_thread_num();
        unsigned T = omp_get_num_threads();
        for (size_t index = t; index < n; index += T)
        {
            partical_sums[t] += vector[index];
        }
    }
    for (size_t i = 1; i < omp_get_num_procs(); ++i)
    {
        partical_sums[0] += partical_sums[i];
    }
    sum = partical_sums[0];
    delete partical_sums;
    return sum / n;
}

double average_omp_with_array_optimized(const double* vector, size_t n)
{
    double sum = 0.0;
    double* partial_sums;
    unsigned T;
#pragma omp parallel
    {
        unsigned t = omp_get_thread_num();
        #pragma omp single 
        {
            T = omp_get_num_threads();
            partial_sums = (double*)malloc(T * sizeof(vector[0]));
        }
        partial_sums[t] = 0.0;
        for (size_t index = t; index < n; index += T)
        {
            partial_sums[t] += vector[index];
        }
    }
    for (size_t i = 1; i < omp_get_num_procs(); ++i)
    {
        partial_sums[0] += partial_sums[i];
    }
    sum = partial_sums[0];
    delete partial_sums;
    return sum / n;
}

double average_omp_with_struct(const double* vector, size_t n)
{
    double sum = 0.0;
    unsigned T;
    partial_sum_t* partial_sums;
#pragma omp parallel
    {
        unsigned t = omp_get_thread_num();
    #pragma omp single
        {
            T = omp_get_num_threads();
            partial_sums = (partial_sum_t*)malloc(T * sizeof(partial_sum_t));
        }
        partial_sums[t].value = 0;
        for (size_t index = t; index < n; index += T)
        {
            partial_sums[t].value += vector[index];
        }
    }
    for (size_t i = 1; i < T; ++i)
    {
        partial_sums[0].value += partial_sums[i].value;
    }
    sum = partial_sums[0].value;
    delete partial_sums;
    return sum / n;
}

double average_omp_CPP(const double* vector, size_t n)
{
    double sum = 0.0;
    double res = 0.0;
#pragma omp parallel
    {
        int t = omp_get_thread_num();
        int T = omp_get_num_threads();
        for (size_t i = t; i < n; i += T)
        {
        #pragma omp critical
            {
                res += vector[i];
            }
        }
    }
    return res / n;
}


double average_omp_mtx(const double* vector, size_t n)
{
    double sum = 0.0;
    double res = 0.0;
    #pragma omp parallel
    {
        int t = omp_get_thread_num();
        int T = omp_get_num_threads();
        double partial_result = 0.0;
        for (size_t i = t; i < n; i += T)
        {
            partial_result += vector[i];
        }
        #pragma omp critical
        {
            res += partial_result;
        }
    }
    return res / n;
}*/