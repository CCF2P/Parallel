#pragma once

#include <vector>
#include <random>
#include <ctime>
#include <iostream>
#include <thread>
#include <mutex>
#include <chrono>

class parallel
{
private:
	std::uniform_real_distribution<double> dist;
	std::mt19937 generator;

	unsigned int T;

public:
	parallel();
	~parallel();

	//////////////////////////////////////////////////////
	/////////////// Methods for benchmark ////////////////
	//////////////////////////////////////////////////////
	std::chrono::time_point<std::chrono::steady_clock> get_time_point();

	//////////////////////////////////////////////////////
	////////////// Array Generation Methods //////////////
	//////////////////////////////////////////////////////
	std::vector<float> generate1D_float(int size, float probability);
	std::vector<double> generate1D_double(int size, float probability);
	std::vector<int> generate1D_int(int size, float probability);

	std::vector<std::vector<float>> generate2D(int size, int probability);
	std::vector<std::vector<float>> generate2D(int width, int height, int probability);

	//////////////////////////////////////////////////////
	////////////////// Simple functions //////////////////
	//////////////////////////////////////////////////////
	float sum(std::vector<float>* data);
	double sum(std::vector<double>* data);
	int sum(std::vector<int>* data);

	//////////////////////////////////////////////////////
	///////////////// Sorting algorithms /////////////////
	//////////////////////////////////////////////////////
	// void odd_even_sort(std::vector<double>* A, int n);

	//////////////////////////////////////////////////////
	// Algorithms for searching a string in a substring //
	//////////////////////////////////////////////////////

	//////////////////////////////////////////////////////
	/////////////////// Other methods ////////////////////
	//////////////////////////////////////////////////////
	template<typename T, typename A>
	void print_array1D(std::vector<T, A>* arr)
	{
		for (int i = 0; i < arr->size(); ++i)
		{
			std::cout << (*arr)[i] << " ";
		}
		std::cout << "\n";
	}

	void print_result(std::chrono::time_point<std::chrono::steady_clock>* t1,
					  std::chrono::time_point<std::chrono::steady_clock>* t2)
	{
		std::chrono::duration<float> d(*t2 - *t1);
		std::cout << "time = "
			<< std::chrono::duration_cast<std::chrono::nanoseconds>(d).count()
			<< " nanoseconds\n";
		std::cout << "time = "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(d).count()
			<< " milliseconds\n";
	}
};