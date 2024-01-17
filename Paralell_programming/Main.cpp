#include <concepts>
#include <condition_variable>
#include <io.h>
#include <iostream>
#include <memory>
#include <mutex>
#include <omp.h>
#include <queue>
#include <ranges> // algorithm - в c++ < 20
#include <stdint.h>
#include <thread>
#include <vector>
#include <functional>
#include <locale>
#include <stdio.h>
#include <string>

#include "m.h"
#include "thread_misc.h"


constexpr unsigned A = 22695477;
constexpr unsigned B = 1;
struct lc_t // linear congruental generator
{
private:
    uint32_t A, B;

public:
    lc_t(uint32_t a = 1, uint32_t b = 0) : A(a), B(b) {}

    lc_t& operator*=(const lc_t x)
    {
        if (A == 1 && B == 0)
        {
            A = x.A;
            B = x.B;
        }
        else
        {
            A *= x.A;
            B += A * x.B;
        }
        return *this;
    }

    auto operator () (uint32_t seed, uint32_t min_val, uint32_t max_val) const
    {
        if (max_val - min_val + 1 == 0)
            return (A * seed + B) % (max_val - min_val) + min_val;
        else
            return A * seed + B;
    }
};


struct partial_sum_t
{
    alignas(64) double value;
};


struct profiling_results_t
{
    double result;
    double time, speedup, efficiency;
    unsigned T;
};


double average_reduce(const uint32_t* V, size_t n)
{
    double res = 0.0;
#pragma omp parallel for reduction(+:res)
    for (int i = 0; i < n; ++i)
    {
        res += V[i];
    }
    return res / n;
}


double average_omp_with_struct(const uint32_t* vector, size_t n)
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
    free(partial_sums);
    return sum / n;
}

double average_cpp_mtx(const uint32_t* v, size_t N)
{
    unsigned T = std::thread::hardware_concurrency();
    std::vector<std::thread> workers;
    std::mutex mtx;

    double res = 0.0;

    auto worker_proc = [&mtx, T, v, N, &res](unsigned t) 
    {
        double partial_result = 0.0;
        for (std::size_t i = t; i < N; i += T)
            partial_result += v[i];

        // mtx.lock();
        std::scoped_lock l{ mtx }; // lock_guard - C++ 11, scopred_lock - C++ 17
        res += partial_result;     // unique_lock
                                 // mtx.unlock();
    };

    for (unsigned t = 0; t < T; ++t)
        workers.emplace_back(worker_proc, t);

    for (auto& w : workers)
        w.join();

    return res / N;
}

double average_cpp_align(const double* v, size_t N)
{
    unsigned T = std::thread::hardware_concurrency();
    std::vector<std::thread> workers;

    double res = 0.0;
    partial_sum_t* partial_sums = (partial_sum_t*)malloc(T * sizeof(partial_sum_t));

    auto worker_proc = [T, v, N, &res, &partial_sums](unsigned t) 
    {
        partial_sums[t].value = 0;
        for (std::size_t i = t; i < N; i += T)
            partial_sums[t].value += v[i];
    };

    for (unsigned t = 0; t < T; ++t)
        workers.emplace_back(worker_proc, t);

    for (auto& w : workers)
        w.join();

    for (size_t i = 1; i < T; ++i)
        partial_sums[0].value += partial_sums[i].value;

    res = partial_sums[0].value;
    free(partial_sums);
    return res / N;
}


double average_cpp_partial_sums_align_local_cach_reduce(const uint32_t* v, size_t N)
{
    unsigned T = std::thread::hardware_concurrency();
    std::vector<std::thread> workers(T - 1);
    std::mutex mtx;

    double res = 0.0;
    std::vector<double> partial_results(T);
    barrier bar(T);

    auto worker_proc = [&mtx, T, v, N, &partial_results, &bar](unsigned t) 
    {
        size_t b = N % T, e = N / T;
        if (t < b)
            b = t * ++e;
        else
            b += t * e;

        e += b;

        double partial_result = 0.0;
        for (size_t i = b; i < e; ++i)
            partial_result += v[i];

        // mtx.lock();
        // std::scoped_lock l { mtx };  // lock_guard - C++ 11, scoped_lock - C++
        // 17
        partial_results[t] = partial_result; // unique_lock
        // mtx.unlock();

        // ====== Reduction begin ======
        for (size_t step = 1, next = 2; step < T; step = next, next += next) {
            bar.arrive_and_wait();
            if ((t & (next - 1)) == 0 && t + step < T) // t % next && t + step < T
            {
                partial_results[t] += partial_results[t + step];
            }
        }
    };

    for (size_t t = 0; t < T - 1; ++t)
        workers[t] = std::thread(worker_proc, t + 1);

    worker_proc(0);
    for (auto& w : workers)
        w.join();

    return partial_results[0] / N;
}


double average_cpp_mtx_align_local_cach(const uint32_t* v, size_t N)
{
    unsigned T = std::thread::hardware_concurrency();
    std::vector<std::thread> workers(T - 1);
    std::mutex mtx;

    double res = 0.0;

    auto worker_proc = [&mtx, T, v, N, &res](unsigned t) 
    {
        size_t b = N % T, e = N / T;
        if (t < b)
            b = t * ++e;
        else
            b += t * e;

        e += b;

        double partial_result = 0.0;
        for (size_t i = b; i < e; ++i)
            partial_result += v[i];

        // mtx.lock();
        std::scoped_lock l { mtx }; // lock_guard - C++ 11, scoped_lock - C++ 17
        res += partial_result;     // unique_lock
                                 // mtx.unlock();
    };

    for (size_t t = 0; t < T - 1; ++t)
        workers[t] = std::thread(worker_proc, t + 1);

    worker_proc(0);
    for (auto& w : workers)
        w.join();

    return res / N;
}


double randomize_vector(uint32_t* V, size_t n, uint32_t seed,
                        uint32_t min_val = 0, uint32_t max_val = UINT32_MAX)
{
    if (min_val > max_val)
        exit(__LINE__);

    double res = 0.0;

    lc_t g = lc_t(A, B);
    for (int i = 1; i < n; ++i)
    {
        g *= g;
        V[i] = g(seed, min_val, max_val);
        res += V[i];
    }

    return res / n;
}


double randomize_vector_thread(
    uint32_t* V,
    size_t n,
    uint32_t seed,
    uint32_t min_val = 0,
    uint32_t max_val = UINT32_MAX)
{
    if (min_val > max_val)
        exit(__LINE__);

    double res = 0.0;
    unsigned T = omp_get_num_threads();
    std::vector<std::thread> workers;
    std::mutex mtx;

    auto worker_proc = [V, n, seed, T, min_val, max_val, &res, &mtx](unsigned t)
    {
        double part_res = 0.0;

        size_t b = n % T, e = n / T;
        if (t < b)
            b = t * ++e;
        else
            b += t * e;
        e += b;

        auto g = fast_pow(lc_t(A, B), fast_pow(2u, b + 1));
        for (size_t i = b; i < e; ++i)
        {
            g *= g;
            V[i] = g(seed, min_val, max_val);
            part_res += V[i];
        }

        std::scoped_lock l{ mtx };
        res += part_res;
    };

    for (size_t t = 1; t < T; ++t)
        workers.emplace_back(worker_proc, t);

    worker_proc(0);
    for (auto& w : workers)
        w.join();

    return res / n;
}

double randomize_vector(std::vector<uint32_t>& v, uint32_t seed,
                        uint32_t max_val = UINT32_MAX, uint32_t min_val = 0)
{
    return randomize_vector(v.data(), v.size(), seed, min_val, max_val);
}


void measure_time(double* f(const double*, size_t), size_t N, std::unique_ptr<double[]>& arr, std::string msg)
{
    double t1 = omp_get_wtime();
    double v = *f(arr.get(), N);
    double t2 = omp_get_wtime();
    std::cout << msg << t2 - t1 << std::endl;
    std::cout << "Result: " << v << std::endl << std::endl;
}

template <class F>
auto measure_time_chrono(F f, size_t N, std::unique_ptr<double[]>& arr)
{
    using namespace std::chrono;
    auto t1 = steady_clock::now();
    f(arr.get(), N);
    auto t2 = steady_clock::now();
    return duration_cast<milliseconds>(t2 - t1).count();
}


void FirstLab(size_t N)
{
    std::vector<uint32_t> buf(N);
    randomize_vector(buf, 20020922);

    std::vector<measure_func> functions_for_measure{
        measure_func("average_reduce", average_reduce),
        measure_func("average_cpp_mtx_align_local_cach", average_cpp_mtx_align_local_cach),
        measure_func("average_cpp_partial_sums_align_local_cach_reduce", average_cpp_partial_sums_align_local_cach_reduce)
    };

    if (_isatty(_fileno(stdout)))
    {
// Код ниже для вывода в консоль
        for (auto& mf : functions_for_measure)
        {
            auto exp_res = run_experiement_cpp(mf.func, N, buf.data());
            std::cout << "Function: " << mf.name << '\n';
            std::cout << "T\tResult\t\t\tTime\t\tSpeedup\t\t\tEfficiency" << '\n';
            for (auto& ev : exp_res)
            {
                std::cout << ev.T << "\t";
                std::cout << ev.result << "\t\t";
                std::cout << ev.time << "\t\t";
                std::cout << ev.speedup << "\t\t\t";
                std::cout << ev.efficiency << '\n';
            }
        }
    }
    else
    {
        std::cout.imbue(std::locale(""));
        // Код ниже для перенаправленного вывода
        std::cout << "Method;T;Result;Time;Speedup;Efficiency\n";
        for (auto& mf : functions_for_measure)
        {
            auto exp_res = run_experiement_cpp(mf.func, N, buf.data());
            for (auto& ev : exp_res)
            {
                std::cout << mf.name << ";";
                std::cout << ev.T << ";";
                std::cout << ev.result << ";";
                std::cout << ev.time << ";";
                std::cout << ev.speedup << ";";
                std::cout << ev.efficiency << "\n";
            }
        }
    }
}


void SecondLab(size_t N)
{
    std::vector<uint32_t> arr(N);

    std::size_t T_max = get_num_threads();
    std::vector<profiling_result_t> profiling_results(T_max);

    for (unsigned T = 1; T <= T_max; ++T)
    {
        set_num_threads(T);

        profiling_results[T - 1].T = T;

        auto t1 = std::chrono::steady_clock::now();
        profiling_results[T - 1].result = randomize_vector(arr, 20020922);
        auto t2 = std::chrono::steady_clock::now();

        profiling_results[T - 1].time = duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        profiling_results[T - 1].speedup = profiling_results[0].time / profiling_results[T - 1].time;
        profiling_results[T - 1].efficiency = profiling_results[T - 1].speedup / T;
    }

    // Вывод результатов
    if (_isatty(_fileno(stdout)))
    {
        // Код ниже для вывода в консоль

        std::cout << "Randomize Vectors" << '\n';
        std::cout << "T\tResult\t\t\tTime\t\tSpeedup\t\t\tEfficiency" << '\n';

        for (auto& pr : profiling_results)
        {
            std::cout << pr.T << "\t";
            std::cout << pr.result << "\t\t";
            std::cout << pr.time << "\t\t";
            std::cout << pr.speedup << "\t\t\t";
            std::cout << pr.efficiency << '\n';
        }
    }
    else
    {
        // Код ниже для перенаправленного вывода
        std::cout.imbue(std::locale(""));
        std::cout << "Method;T;Result;Time;Speedup;Efficiency\n";

        for (auto& pr : profiling_results)
        {
            std::cout << "Randomize Vectors;";
            std::cout << pr.T << ";";
            std::cout << pr.result << ";";
            std::cout << pr.time << ";";
            std::cout << pr.speedup << ";";
            std::cout << pr.efficiency << "\n";
        }
    }
}


int main()
{
    const size_t N1 = 1u << 26;
    const size_t N2 = 1u << 26;
    FirstLab(N1);
    SecondLab(N2);
    return 0;
}