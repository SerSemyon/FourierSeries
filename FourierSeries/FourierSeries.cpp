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
std::normal_distribution<double> distribution(0.0, 0.5);

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

void Normalize(std::vector<int>& x) {
    int carry = 0;
    for (int i = x.size() - 1; i >= 0; --i) {
        x[i] += carry;
        carry = x[i] / 10;
        x[i] %= 10;
    }
}

std::vector<int> Multiply(const std::vector<int>& a, const std::vector<int>& b) {
    std::vector < std::complex < double>> fa(a.rbegin(), a.rend()), fb(b.rbegin(), b.rend());
    size_t n = 1;
    for (; n <= std::max(a.size(), b.size()); n <<= 1) {}
    n <<= 1;

    fa.resize(n), fb.resize(n);

    std::vector<std::complex<double>> resA = fa;
    std::vector<std::complex<double>> resB = fb;

    for (int i = 0; i < n; i++)
    {
        resA[i] = FFT(fa, n, i);
        resB[i] = FFT(fb, n, i);

    }
    for (size_t i = 0; i < n; ++i) {
        resA[i] *= resB[i];
    }
    std::vector<std::complex<double>> out = resA;
    for (int i = 0; i < n; i++)
    {
        out[i] = IFFT(resA, n, i);
    }

    std::vector<int> result(n);
    for (size_t i = 0; i < n; ++i) {
        result[i] = std::round(out[i].real());
    }

    Normalize(result);
    std::reverse(result.begin(), result.end());
    size_t shift = 0;
    for (size_t i = 0; i < result.size(); ++i) {
        bool is_null = result[i] == 0;
        shift += (is_null) ? 1 : 0;
        if (!is_null) {
            break;
        }
    }
    result = { result.begin() + shift, result.end() };
    return result;
}

void MultipleBigInteger()
{
    std::cout << "Multiply big integer" << std::endl;
    unsigned int a_size = 4;
    std::vector<int> a(a_size);
    for (int i = 0; i < a_size; i++)
    {
        a[i] = rand()%10;
        std::cout << a[i];
    }
    std::cout << std::endl;
    unsigned int b_size = 4;
    std::vector<int> b(b_size);
    for (int i = 0; i < b_size; i++)
    {
        b[i] = rand() % 10;
        std::cout << b[i];
    }
    std::cout << std::endl;
    auto result = Multiply(a, b);
    for (int i = 0; i < result.size(); i++)
    {
        std::cout << result[i];
    }
    std::cout << std::endl;
}

double hCore(double omega, double t)
{
    return sin(omega * M_PI * t) / (M_PI * t);
}

template<typename T>
std::vector<T> sincFilter(std::vector<T> x, double omega)
{
    std::vector<T> y = new vector<T>(x.size);
    for (int n = 0; n < x.size; n++)
    {
        y[n] = 0;
        for (int k = 0; k < n; k++)
        {
            y[n] += hCore(omega, k) * x[n - k];
        }
    }
    return y;
}

template<typename T>
std::vector<T> firstCriteria(const std::vector<T>& X, double noiseAmplitude) {
    std::vector<T> out(X.begin(), X.end());
    std::vector<double> amplitude(X.size());
    double maxXk = 0;
    for (size_t k = 0; k < out.size(); k++) {
        amplitude[k] = abs(X[k]);
        //double absXk = abs(out[k]);
        if (amplitude[k] > maxXk)
            maxXk = amplitude[k];
    }
    double epsilon = noiseAmplitude / maxXk + 0.1;
    for (size_t k = 0; k < out.size(); k++) {
        if (amplitude[k] / maxXk <= epsilon)
            out[k] = 0;
    }
    return out;
}

template<typename T>
std::vector<T> secondCriteria(const std::vector<T>& X) {
    std::vector<T> out(X.begin(), X.end());
    std::vector<double> amplitude(X.size());
    double maxXk = pow(abs(out[0]), 2);
    double minXk = maxXk;
    double sum = 0;
    for (size_t k = 0; k < out.size(); k++) {
        amplitude[k] = pow(abs(out[k]), 2);
        if (amplitude[k] > maxXk)
            maxXk = amplitude[k];
        if (amplitude[k] < minXk)
            minXk = amplitude[k];
        sum += amplitude[k];
    }
    sum /= out.size();
    double epsilon = minXk / sum + 0.1;
    for (size_t k = 0; k < out.size(); k++) {
        if (amplitude[k] / maxXk <= epsilon)
            out[k] = 0;
    }
    return out;
}

void FiltrTask1() 
{
    std::cout << "First criteria" << std::endl;
    unsigned int n = 100;
    std::vector<double> x(n);
    double* y = new double[n];
    std::vector<std::complex<double>> complexX(n);
    std::vector<std::complex<double>> complexY(n);
    std::cout << "x = [" << std::endl;
    for (int i = 0; i < n; i++)
    {
        x[i] = 0.01 * i;
        y[i] = X_2(x[i]);
        complexX[i] = { y[i],0 };
        std::cout << y[i] << ",";
    }
    std::cout << "];" << std::endl;
    for (int i = 0; i < n; i++)
    {
        complexY[i] = FFT(complexX, n, i);
    }
    std::cout << std::endl;
    complexY = firstCriteria<std::complex<double>>(complexY, 1.5);
    std::cout << "y = [";
    for (int i = 0; i < n; i++)
    {
        std::cout << IFFT(complexY, n, i).real() << ",";
    }
    std::cout << "];" << std::endl;
}

void FiltrTask2()
{
    std::cout << "Second criteria" << std::endl;
    unsigned int n = 100;
    std::vector<double> x(n);
    double* y = new double[n];
    std::vector<std::complex<double>> complexX(n);
    std::vector<std::complex<double>> complexY(n);
    std::cout << "x = [" << std::endl;
    for (int i = 0; i < n; i++)
    {
        x[i] = 0.03 * i;
        y[i] = X_2(x[i]);
        complexX[i] = { y[i],0 };
        std::cout << y[i] << ",";
    }
    std::cout << "];" << std::endl;
    for (int i = 0; i < n; i++)
    {
        complexY[i] = FFT(complexX, n, i);
    }
    std::cout << std::endl;
    complexY = secondCriteria<std::complex<double>>(complexY);
    std::cout << "y = [";
    for (int i = 0; i < n; i++)
    {
        std::cout << IFFT(complexY, n, i).real() << ",";
    }
    std::cout << "];" << std::endl;
}

int main()
{

    //Test();
    //Task4();
    ////Task5();
    FiltrTask2();

}
