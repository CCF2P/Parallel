#pragma once
#include <iostream>
#include <vector>

std::vector<int> prefix_function(const std::string& s)
{
	std::vector<int> pi(s.length(), 0);
	int j = 0;
	for (int i = 1; i < s.length(); i++)
	{
		j = pi[i - 1];  //текущая длина префикса, который мы хотим продолжить
		//гарантируется, что s[0..j-1] = s[i-j..i-1].

		while (j > 0 && s[i] != s[j]) //пока мы не можем продолжить текущий префикс
		{
			j = pi[j - 1];  //уменьшаем его длину до следующей возможной
		}

		//Теперь j - максимальная длина префикса, который мы можем продолжить,
		//или 0, если такового не существует.

		if (s[i] == s[j])
		{
			pi[i] = j + 1;
		}
		else
		{    //такое может произойти только при j = 0
			pi[i] = j;
		}
	}

	return pi;
}

std::vector<int> KMP(const std::string& text, const std::string& pattern)
{
	int n = text.length();
	int m = pattern.length();
	std::vector<int> pi = prefix_function(text);
	std::vector<int> indexes;

	int i = 0;
	int j = 0;
	while (i < n)
	{
		if (pattern[j] == text[i])
		{
			i++;
			j++;
		}

		if (j == m)
		{
			indexes.emplace_back(i - m);
			j = pi[j - 1];
		}
		else if (i < n && pattern[j] != text[i])
		{  // If there is a mismatch
			if (j == 0)          // if j becomes 0 then simply increment the index i
				i++;
			else
				j = pi[j - 1];  //Update j as Lps of last matched character
		}
	}
	return indexes;
}

static bool compare(double a, double b)
{
	if (a > b)
		return 1;
	return 0;
}

static void exchange(std::vector<double>* A, int a_idx, int b_idx)
{
	double tmp = (*A)[a_idx];
	(*A)[a_idx] = (*A)[b_idx];
	(*A)[b_idx] = tmp;
}

static void OddEvenSort(std::vector<double>* A, int n)
{
	for (int i = 0; i < n; i++)
	{
		if (i % 2 == 1)
		{
			for (int j = 0; j < n / 2 - 1; ++j)
			{
				int t1 = 2 * j + 1;
				int t2 = 2 * j + 2;
				if (compare(A->data()[t1], A->data()[t2]))
					exchange(A, t1, t2);
			}
		}
		else
		{
			for (int j = 0; j < n / 2; ++j)
			{
				int t1 = 2 * j;
				int t2 = 2 * j + 1;
				if (compare(A->data()[t1], A->data()[t2]))
					exchange(A, t1, t2);
			}
		}
	}
}