#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <complex>
#include <vector>
#include <random>
#include <ctime>

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
    std::cout << "input" << std::endl;
    for (int i = 0; i < 4; i++)
    {
        std::cout << x[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "after DFT" << std::endl;
    for (int i = 0; i < 4; i++)
    {
        complex<double> current_point = DFT(x, 4, i);
        res[i] = current_point;
        std::cout << current_point << std::endl;
    }
    std::cout << std::endl;
    std::cout << "after IDFT" << std::endl;
    for (int i = 0; i < 4; i++)
    {
        std::cout << IDFT(res, 4, i) << std::endl;
    }
    std::cout << std::endl;
}

double X_1(double n)
{
    return sin(2 * n) + sin(3 * n) + sin(5 * n);
}

void Task4() 
{
    std::cout << "Task 4" << std::endl;
    unsigned int n = 100;
    std::vector<double> x(n);
    double* y = new double[n];
    std::vector<std::complex<double>> complexX(n);
    std::vector<std::complex<double>> complexY(n);
    std::cout << "input" << std::endl;
    for (int i = 0; i < n; i++)
    {
        x[i] = 0.01 * i;
        y[i] = X_1(x[i]);
        complexX[i] = { y[i],0 };
        std::cout << x[i] << " " << complexX[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "after DFT" << std::endl;
    for (int i = 0; i < n; i++)
    {
        complexY[i] = DFT(complexX, n, i);
        std::cout << x[i] << " " << complexY[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "after IDFT" << std::endl;
    for (int i = 0; i < n; i++)
    {
        std::cout << x[i] << " " << IDFT(complexY, n, i) << std::endl;
    }
    std::cout << std::endl;
}

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0, 1.0);

double X_2(double x)
{
    double x1 = X_1(x);
    double res = x1 + distribution(generator);
    return res;
}

void Task5() {
    std::cout << "Task 5" << std::endl;
    unsigned int n = 100;
    std::vector<double> x(n);
    double* y = new double[n];
    std::vector<std::complex<double>> complexX(n);
    std::vector<std::complex<double>> complexY(n);
    std::cout << "input" << std::endl;
    for (int i = 0; i < n; i++)
    {
        x[i] = 0.01 * i;
        y[i] = X_2(x[i]);
        complexX[i] = { y[i],0 };
        std::cout << x[i] << " " << complexX[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "after DFT" << std::endl;
    for (int i = 0; i < n; i++)
    {
        complexY[i] = DFT(complexX, n, i);
        std::cout << x[i] << " " << complexY[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "after IDFT" << std::endl;
    for (int i = 0; i < n; i++)
    {
        std::cout << x[i] << " " << IDFT(complexY, n, i) << std::endl;
    }
    std::cout << std::endl;
}

complex<double> FFT(std::vector<std::complex<double>> x, unsigned int N, double k)
{
    complex<double> W_N = E_ix(2 * M_PI / N);
    complex<double> X_1 = 0;
    complex<double> X_2 = 0;
    for (int m = 0; m < N / 2; m++)
    {
        complex<double> W_Nnk = 1;
        for (int i = 0; i < 2 * m * k; i++)
        {
            W_Nnk *= W_N;
        }
        X_1 += x[2 * m] * W_Nnk;
        X_2 += x[2 * m + 1] * W_Nnk;
    }
    complex<double> W_Nk = 1;
    for (int i = 0; i < k; i++)
    {
        W_Nk *= W_N;
    }
    return X_1 + X_2 * W_Nk;
}

complex<double> IFFT(std::vector<std::complex<double>> x, unsigned int N, double k)
{
    complex<double> W_N = E_ix(2 * M_PI / N);
    complex<double> X_1 = 0;
    complex<double> X_2 = 0;
    for (int m = 0; m < N / 2; m++)
    {
        complex<double> W_Nnk = 1;
        for (int i = 0; i < 2 * m * k; i++)
        {
            W_Nnk /= W_N;
        }
        X_1 += x[2 * m] * W_Nnk;
        X_2 += x[2 * m + 1] * W_Nnk;
    }
    X_1/=N;
    X_2/=N;
    complex<double> W_Nk = 1;
    for (int i = 0; i < k; i++)
    {
        W_Nk /= W_N;
    }
    return X_1 + X_2 * W_Nk;
}

void Task1() {
    std::cout << "Task 1" << std::endl;
    unsigned int n = 10;
    std::vector<double> x(n);
    double* y = new double[n];
    std::vector<std::complex<double>> complexX(n);
    std::vector<std::complex<double>> complexY(n);
    std::cout << "input" << std::endl;
    for (int i = 0; i < n; i++)
    {
        x[i] = 0.01 * i;
        y[i] = X_2(x[i]);
        complexX[i] = { y[i],0 };
        std::cout << x[i] << " " << complexX[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "after DFT" << std::endl;
    for (int i = 0; i < n; i++)
    {
        complexY[i] = FFT(complexX, n, i);
        std::cout << x[i] << " " << complexY[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "after IDFT" << std::endl;
    for (int i = 0; i < n; i++)
    {
        std::cout << x[i] << IFFT(complexY, n, i) << std::endl;
    }
    std::cout << std::endl;
}

unsigned int FindExecutionTime(void method())
{
    unsigned int start_time = clock(); // начальное время
    method();
    unsigned int end_time = clock(); // конечное время
    unsigned int search_time = end_time - start_time;
    return search_time;
}

int main()
{
    //Test();
    //Task4();
    //Task5();
    //std::cout << FindExecutionTime(Task5) << std::endl;
    std::cout << FindExecutionTime(Task1) << std::endl;
}
