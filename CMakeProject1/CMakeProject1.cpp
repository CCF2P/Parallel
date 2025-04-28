// CMakeProject1.cpp: определяет точку входа для приложения.
//

#include "CMakeProject1.h"

int main()
{
	parallel p;

	auto t1 = p.get_time_point();
	std::vector<int> vec = p.generate1D_int(2 << 10, 0.9);
	auto t2 = p.get_time_point();
	p.print_result(&t1, &t2);
	p.print_array1D(&vec);

	return 0;
}
