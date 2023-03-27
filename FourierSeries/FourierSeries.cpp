#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <complex>
#include <vector>

using namespace std;

complex<double> E_ix(double x)
{
    return complex<double>(cos(x), -sin(x));
}

complex<double> DFT(std::vector<std::complex<double>> x, unsigned int N, double k)
{
    complex<double> W_N = E_ix(2*M_PI/N);
    complex<double> sum = 0;
    for (int n = 0; n < N; n++)
    {
        complex<double> W_Nnk = 1;
        for (int i = 0; i < n * k; i++)
        {
            W_Nnk *= W_N;
        }
        sum += x[n] * W_Nnk;
    }
    return sum;
}

complex<double> IDFT(std::vector<std::complex<double>> x, unsigned int N, double n)
{
    complex<double> W_N = E_ix(2 * M_PI / N);
    complex<double> sum = 0;
    for (int k = 0; k < N; k++)
    {
        complex<double> W_Nnk = 1;
        for (int i = 0; i < n * k; i++)
        {
            W_Nnk /= W_N;
        }
        sum += x[k] * W_Nnk;
    }
    sum /= N;
    return sum;
}

void Test()
{
    std::vector<std::complex<double>> x = {
        {1, 0},
        {2, 0},
        {3, 0},
        {4, 0}
    };
    std::vector<std::complex<double>> res = {
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0}
    };
    for (int i = 0; i < 4; i++)
    {
        complex<double> current_point = DFT(x, 4, i);
        res[i] = current_point;
        std::cout << current_point << std::endl;
    }
    std::cout << std::endl;
    for (int i = 0; i < 4; i++)
    {
        std::cout << IDFT(res, 4, i) << std::endl;
    }
}

int main()
{
    Test();
}
