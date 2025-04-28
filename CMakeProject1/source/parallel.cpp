#include "../include/parallel.h"

parallel::parallel()
{
	this->dist = std::uniform_real_distribution<double>(0.0, 1.0);
	this->generator = std::mt19937(std::random_device{}());

	this->T = std::thread::hardware_concurrency();
}

parallel::~parallel() {}

std::chrono::time_point<std::chrono::steady_clock> parallel::get_time_point()
{
	std::chrono::time_point<std::chrono::steady_clock> t = std::chrono::steady_clock::now();
	return t;
}

//////////////////////////////////////////////////////
////////////// Array Generation Methods //////////////
//////////////////////////////////////////////////////
std::vector<float> parallel::generate1D_float(int size, float probability)
{
	std::vector<float> vec(size, 0.0);
	std::vector<std::thread> workers;

	auto worker_proc = [&vec, probability, size](
		unsigned int T, unsigned int t,
		std::mt19937 generator, std::uniform_real_distribution<double> dist
	)
	{
		size_t b = size % T, e = size / T;
		if (t < b)
			b = t * ++e;
		else
			b += t * e;
		e += b;

		float tmp = 0.0f;
		for (size_t i = b; i < e; ++i)
		{
			tmp = dist(generator);
			if (tmp < probability)
			{
				vec[i] = (float)(tmp * 10);
			}
		}
	};

	// Запускаем потоки
	for (unsigned int t = 0; t < this->T; ++t)
	{
		workers.emplace_back(std::thread(worker_proc, this->T, t,
							 this->generator, this->dist));
	}

	// Дожидаемся их работы
	for (unsigned int t = 0; t < this->T; ++t)
	{
		workers[t].join();
	}

	return vec;
}

std::vector<double> parallel::generate1D_double(int size, float probability)
{
	std::vector<double> vec(size, 0);
	std::vector<std::thread> workers;

	auto worker_proc = [&vec, probability, size](
		unsigned int T, unsigned int t,
		std::mt19937 generator, std::uniform_real_distribution<double> dist
	)
	{
		size_t b = size % T, e = size / T;
		if (t < b)
			b = t * ++e;
		else
			b += t * e;
		e += b;

		float tmp = 0.0f;
		for (size_t i = b; i < e; ++i)
		{
			tmp = dist(generator);
			if (tmp < probability)
			{
				vec[i] = (double)(tmp * 10);
			}
		}
	};

	// Запускаем потоки
	for (unsigned int t = 0; t < this->T; ++t)
	{
		workers.emplace_back(std::thread(worker_proc, this->T, t,
							 this->generator, this->dist));
	}

	// Дожидаемся их работы
	for (unsigned int t = 0; t < this->T; ++t)
	{
		workers[t].join();
	}

	return vec;
}

std::vector<int> parallel::generate1D_int(int size, float probability)
{
	std::vector<int> vec(size, 0);
	std::vector<std::thread> workers;

	auto worker_proc = [&vec, probability, size](
		unsigned int T, unsigned int t,
		std::mt19937 generator, std::uniform_real_distribution<double> dist
	)
	{
		size_t b = size % T, e = size / T;
		if (t < b)
			b = t * ++e;
		else
			b += t * e;
		e += b;

		float tmp = 0.0f;
		for (size_t i = b; i < e; ++i)
		{
			tmp = dist(generator);
			if (tmp < probability)
			{
				vec[i] = (int)(tmp * 10);
			}
		}
	};

	// Запускаем потоки
	for (unsigned int t = 0; t < this->T; ++t)
	{
		workers.emplace_back(std::thread(worker_proc, this->T, t,
							 this->generator, this->dist));
	}

	// Дожидаемся их работы
	for (unsigned int t = 0; t < this->T; ++t)
	{
		workers[t].join();
	}

	return vec;
}

std::vector<std::vector<float>> parallel::generate2D(int size, int probability)
{
	std::vector<std::vector<float>> mtrx(size, std::vector<float>(size, 0));
	float tmp = 0.0f;

	for (unsigned i = 0; i < size; ++i)
	{
		for (unsigned j = 0; j < size; ++j)
		{
			tmp = this->dist(this->generator);
			if (tmp < probability)
				mtrx[i][j] = tmp * 10;
		}
	}

	return mtrx;
}

std::vector<std::vector<float>> parallel::generate2D(int width, int height, int probability)
{
	std::vector<std::vector<float>> mtrx(width, std::vector<float>(height, 0));
	float tmp = 0.0f;

	for (unsigned i = 0; i < width; ++i)
	{
		for (unsigned j = 0; j < height; ++j)
		{
			tmp = this->dist(this->generator);
			if (tmp < probability)
				mtrx[i][j] = tmp * 10;
		}
	}

	return mtrx;
}

//////////////////////////////////////////////////////
////////////////// Simple functions //////////////////
//////////////////////////////////////////////////////

static void work_fsum(unsigned int T, std::vector<float>* result,
			   unsigned int t, int N, std::vector<float>* data)
{
	float partial_sum = 0.0f;

	for (unsigned i = t; i < N; i += T)
	{
		partial_sum += (*data)[i];
	}

	(*result)[t] = partial_sum;
}

float parallel::sum(std::vector<float>* data)
{
	const int N = data->size();
	std::vector<float> partial_sum(this->T, 0);
	std::vector<std::thread> workers;

	// Запускаем потоки
	for (unsigned int t = 0; t < this->T; ++t)
	{
		workers.emplace_back(std::thread(work_fsum, this->T, &partial_sum,
							 t, N, data));
	}

	// Дожидаемся их работы
	for (unsigned int t = 0; t < this->T; ++t)
	{
		workers[t].join();
	}

	float sum = 0.0f;
	for (unsigned int i = 0; i < this->T; ++i)
	{
		sum += partial_sum[i];
	}

	return sum;
}

static void work_dsum(unsigned int T, std::vector<double>* result,
					  unsigned int t, int N, std::vector<double>* data)
{
	int partial_sum = 0;

	for (unsigned i = t; i < N; i += T)
	{
		partial_sum += (*data)[i];
	}

	(*result)[t] = partial_sum;
}

double parallel::sum(std::vector<double>* data)
{
	const int N = data->size();
	std::vector<double> partial_sum(this->T, 0.0);
	std::vector<std::thread> workers;

	// Запускаем потоки
	for (unsigned int t = 0; t < this->T; ++t)
	{
		workers.emplace_back(std::thread(work_dsum, this->T, &partial_sum,
							 t, N, data));
	}

	// Дожидаемся их работы
	for (unsigned int t = 0; t < this->T; ++t)
	{
		workers[t].join();
	}

	double sum = 0.0;
	for (unsigned int i = 0; i < this->T; ++i)
	{
		sum += partial_sum[i];
	}

	return sum;
}

static void work_sum(unsigned int T, std::vector<int>* result,
			   unsigned int t, int N, std::vector<int>* data)
{
	int partial_sum = 0;

	for (unsigned i = t; i < N; i += T)
	{
		partial_sum += (*data)[i];
	}

	(*result)[t] = partial_sum;
}

int parallel::sum(std::vector<int>* data)
{
	const int N = data->size();
	std::vector<int> partial_sum(this->T, 0);
	std::vector<std::thread> workers;

	// Запускаем потоки
	for (unsigned int t = 0; t < this->T; ++t)
	{
		workers.emplace_back(std::thread(work_sum, this->T, &partial_sum,
							 t, N, data));
	}

	// Дожидаемся их работы
	for (unsigned int t = 0; t < this->T; ++t)
	{
		workers[t].join();
	}

	int sum = 0;
	for (unsigned int i = 0; i < this->T; ++i)
	{
		sum += partial_sum[i];
	}

	return sum;
}

//////////////////////////////////////////////////////
///////////////// Sorting algorithms /////////////////
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
// Algorithms for searching a string in a substring //
//////////////////////////////////////////////////////
/*
static void work_KMP(int* pattern, int* target, int f[], int c[], int n, int m)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int i = n * index;
	int j = n * (index + 2) - 1;
	if (i > m)
		return;
	if (j > m)
		j = m;
	int k = 0;
	while (i < j)
	{
		if (k == -1)
		{
			i++;
			k = 0;
		}
		else if (target[i] == pattern[k])
		{
			i++;
			k++;
			if (k == n)
			{
				c[i - n] = i - n;
				i = i - k + 1;
			}
		}
		else
		{
			k = f[k];
		}
	}
}

std::vector<int> parallel::KMP()
{
	std::vector<int> indexes;



	return indexes;
}*/