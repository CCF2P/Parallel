#include <concepts>
#include <condition_variable>
#include <io.h>
#include <iostream>
#include <memory>
#include <mutex>
#include <omp.h>
#include <queue>
#include <stdint.h>
#include <thread>
#include <vector>




/*
double get_omp_time(double (double (*f)(const double*, size_t), double const* V,
size_t N))
{
    auto t1 = omp_get_wtime();
    f(V, N);
    return omp_get_wtime() - t1;
}*/
/*
template<class F> int requires { std::is_invocable<F> }
auto cpp_get_time(F f, const double* V, size_t N)
{

}
*/
/*
template<std::invocable<const double*, size_t> F, const double*, size_t>
auto cpp_get_time(F f, const double* V, size_t N)
{
    using namespace std::chrono;
    auto t1 = steady_clock::now();
    f(V, N);
    auto t2 = steady_clock::now();
    return duration_cast<millisecond>(t2 - t1).count();
}*/

/*
profiling_results_t* run_experiment(double (*f)(const double* V, size_t n),
const double* V, size_t n)
{
    profiling_results_t* res_table =
(profiling_results_t*)malloc(omp_get_num_procs() * sizeof(profiling_results_t));
    for (unsigned T = 1; T <= omp_get_num_procs(); ++T)
    {
        omp_set_num_threads(T);
        auto t1 = omp_get_wtime();
        res_table[T - 1].result = f(V, n);
        auto t2 = omp_get_wtime();

        res_table[T - 1].time = t2 - t1;
        res_table[T - 1].speedup = res_table[0].time / (t2 - t1);
        res_table[T - 1].efficiency = res_table[T - 1].speedup / T;
        res_table[T - 1].T = T;
    }
    return res_table;
}*/



template<class T, std::unsigned_integral U>
auto my_pow(T x, U n)
{
  T r = T(1);
  while (n > 0) {
    if (n & 1)
      r *= x;
    x *= x;
    x >>= 1;
  }
  return r;
}
/*
double randomize_vector(uint32_t* V,
                 size_t n,
                 uint32_t seed,
                 uint32_t min_val = 0,
                 uint32_t max_val = UINT32_MAX)
{
  constexpr unsigned A = 1;
  constexpr unsigned B = 1;
  double res = 0.0;

  for (int i = 0; i < n; ++i) {
    seed = A * seed + B;
    V[i] = min_val + seed % (max_val - min_val + 1);
    res += V[i];
  }

  return res / n;
}*/
