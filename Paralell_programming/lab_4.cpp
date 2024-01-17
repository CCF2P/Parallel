#include <concepts>
#include <condition_variable>
#include <io.h>
#include <iostream>
#include <memory>
#include <mutex>
#include <omp.h>
#include <queue>
#include <ranges> // algorithm - â c++ < 20
#include <stdint.h>
#include <thread>
#include <vector>


template<class T, std::unsigned_integral U>
auto fast_pow(T x, U n) requires requires(T x) { T(1); x *= x; }
{
    T r = T(1);
    while (n > 0) {
        if (n & 1)
            r *= x;
        x *= x;
        n >>= 1;
    }
    return r;
}


/*int main()
{
    int n = 10;
    uint32_t* V = new uint32_t[n];

    double res = randomize_vector_thread(V, n, 0, 1, 10);
    std::cout << "res = " << res << "\n";
    for (int i = 0; i < n; ++i) 
    {
        std::cout << V[i] << " ";
    }

    std::cout << "\n test:\t" << randomize_test() << "\n";

    delete[] V;
    return 0;
}*/