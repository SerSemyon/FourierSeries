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

double hCore(double f0, double t)
{
    return sin(2* M_PI * f0 * t) / (M_PI * t);
}

template<typename T>
std::vector<T> sincFilter(const std::vector<T>& X, double omega)
{
    std::vector<T> out(X.begin(), X.end());
    for (int n = 0; n < X.size(); n++)
    {
        out[n] = 0;
        for (int k = 1; k < n; k++)
        {
            out[n] += hCore(omega, k) * X[n - k];
        }
    }
    return out;
}

void FiltrTask3()
{
    std::cout << "Sinc filtr" << std::endl;
    unsigned int n = 100;
    std::vector<double> x(n);
    double* y = new double[n];
    std::vector<std::complex<double>> complexX(n);
    std::vector<std::complex<double>> complexY(n);
    std::cout << "x = [" << std::endl;
    for (int i = 0; i < n; i++)
    {
        x[i] = 0.03 * i+0.1;
        y[i] = X_1(x[i]);
        complexX[i] = { y[i],0 };
        std::cout << y[i] << ",";
    }
    std::cout << "];" << std::endl;
    for (int i = 0; i < n; i++)
    {
        complexY[i] = FFT(complexX, n, i);
    }
    std::cout << std::endl;
    complexY = sincFilter<std::complex<double>>(complexX, 10);
    std::cout << "y = [";
    for (int i = 0; i < n; i++)
    {
        //std::cout << IFFT(complexY, n, i).real() << ",";
        std::cout << complexY[i].real() << ",";
    }
    std::cout << "];" << std::endl;
}

double thetaBessel(double x, int n)
{
    unsigned int factNpK = 1;
    unsigned int factNmK = 1;
    unsigned int factK = 1;
    for (int i = 0; i < n; i++)
        factNpK *= i;
    factNmK = factNpK;
    double sum = pow(x, n);
    for (int k = 1; k < n; k++)
    {
        factK *= k;
        factNpK *= (n + k);
        factNmK *= (n - k);
        sum += pow(x, n - k) * factNpK / (factNmK * factK * pow(2, k));
    }
    return sum;
}

template<typename T>
std::vector<T> BesselFilter(std::vector<T>& X, double omega0)
{
    std::vector<T> out(X.begin(), X.end());
    int n = 4;
    for (size_t k = 0; k < out.size(); k++) {
        out[k] = X[k] * thetaBessel(0, n) / thetaBessel(abs(X[k]) / omega0, n);
    }
    return out;
}

void FiltrTask4()
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
        //complexY[i] = FFT(complexX, n, i);
    }
    std::cout << std::endl;
    complexY = BesselFilter<std::complex<double>>(complexX, 10);
    std::cout << "y = [";
    for (int i = 0; i < n; i++)
    {
        std::cout << IFFT(complexY, n, i).real() << ",";
    }
    std::cout << "];" << std::endl;
}

template <typename Iterator>
std::vector<std::complex<double>> Fft(Iterator begin, Iterator end, bool invert = false) {
    std::vector<std::complex<double>> x(begin, end);
    if (x.size() != 1 && std::fmod(std::log2(x.size()), 1) != 0) {
        throw std::invalid_argument("Размерность сигнала должна быть степенью двойки ");
    }

    size_t n = x.size();
    if (n == 1)  return x;

    std::vector<std::complex<double>> x1(n / 2), x2(n / 2);
    for (size_t i = 0, j = 0; i < n; i += 2, ++j) {
        x1[j] = x[i];
        x2[j] = x[i + 1];
    }
    x1 = Fft(x1.begin(), x1.end(), invert);
    x2 = Fft(x2.begin(), x2.end(), invert);

    double theta = 2 * M_PI / n * (invert ? -1 : 1);
    std::complex<double> w(1), wn(std::cos(theta), std::sin(theta));
    for (size_t i = 0; i < n / 2; ++i) {
        x[i] = x1[i] + w * x2[i];
        x[i + n / 2] = x1[i] - w * x2[i];
        if (invert) {
            x[i] /= 2, x[i + n / 2] /= 2;
        }
        w *= wn;
    }
    return x;
}

void WaveTask1() {
    std::cout << "WaveTask1" << std::endl;
    size_t size = 64;
    double h = 0.2;
    std::vector<double> x(size);
    for (size_t i = 0; i < size; i++)
    {
        x[i] = h * i;
    }
    std::vector<double> y(size);
    for (size_t i = 0; i < size; i++) {
        y[i] = std::atan(2 * x[i] - 6.5) * x[i] + x[i] * std::sin(x[i] + std::cos(2 * x[i]));
    }
    for (size_t i = 0; i < y.size(); ++i) {
        std::cout << x[i] << ',' << y[i] << '\n';
    }
    std::cout << "\n\n";
    size_t size_window = 32;
    std::vector<std::vector<std::complex<double>>> vec_X;
    vec_X.reserve(std::ceil(size / size_window) + 1);
    for (auto beg = x.begin(), end = x.begin() + size_window; end < x.end(); ++beg, ++end) {

        auto beg_y = y.begin() + std::distance(x.begin(), beg);
        auto end_y = beg_y + size_window;
        vec_X.push_back(Fft(beg_y, end_y));
    }
    for (size_t i = 0; i < vec_X.size(); ++i) {
        for (size_t j = 0; j < vec_X[i].size(); ++j) {
            std::cout << i << ',' << j << ',' << std::abs(vec_X[i][j]) << '\n';
        }
    }
    std::cout << std::endl << std::endl;
}

double Haar(double x) {
    if (x < 1 && x >= 0.5) {
        return -1;
    }
    else if (x > 0 && x < 0.5) {
        return 1;
    }
    return 0;
}

double MHAT(double x) {
    return (1 - x * x) * std::exp(-(x * x) / 2.0);
}

void ExampleWaveletMHAT() {
    std::cout << "ExampleWaveletMHAT" << std::endl;
    size_t size = 64;
    double b0 = 1, a0 = 0.5;
    double h = 0.2;
    std::vector<double> x(size);
    for (size_t i = 0; i < size; i++)
    {
        x[i] = h * i;
    }
    std::vector<double> mhat(size);
    std::vector<double> y(size);

    auto psi = [b0, a0](double x, size_t m, size_t n) {
        return (1.0 / std::sqrt(m)) * MHAT((x - n) / m);
    };

    for (size_t i = 0; i < size; ++i) {
        mhat[i] = MHAT(x[i]);
        y[i] = std::sin(2 * x[i]) + std::sin(3 * x[i]);
    }
    std::vector<std::vector<double>> W(x.size(), std::vector<double>(x.size(), 0));
    for (size_t i = 0; i < mhat.size(); ++i) {
        std::cout << x[i] << ',' << mhat[i] << '\n';
    }
    std::cout << "\n\n";
    for (size_t i = 0; i < y.size(); ++i) {
        std::cout << x[i] << ',' << y[i] << '\n';
    }
    std::cout << "\n\n";
    for (size_t m = 0; m < x.size(); ++m) {
        for (size_t n = 0; n < x.size(); ++n) {
            double scalar_product = 0.0;
            // вычисляем скалярное произведение
            for (size_t i = 0; i < x.size(); ++i) {
                scalar_product += y[i] * psi(i, m, n);
            }
            W[m][n] = scalar_product;
        }
    }
    for (size_t i = 0; i < W.size(); ++i) {
        for (size_t j = 0; j < W[i].size(); ++j) {
            std::cout << i << ',' << j << ',' << W[i][j] << '\n';
        }
    }
    std::cout << "\n\n";
}

void ExampleWaveletHAAR() {
    std::cout << "ExampleWaveletHAAR" << std::endl;
    size_t size = 64;
    double b0 = 1, a0 = 0.5;
    double h = 0.2;
    std::vector<double> x(size);
    for (size_t i = 0; i < size; i++)
    {
        x[i] = h * i;
    }
    std::vector<double> mhat(size);
    std::vector<double> y(size);

    auto psi = [b0, a0](double x, size_t m, size_t n) {
        return (1.0 / std::sqrt(m)) * Haar((x - n) / m);
    };

    for (size_t i = 0; i < size; ++i) {
        mhat[i] = Haar(x[i]);
        y[i] = std::sin(2 * x[i]) + std::sin(3 * x[i]);
    }
    std::vector<std::vector<double>> W(x.size(), std::vector<double>(x.size(), 0));
    for (size_t i = 0; i < mhat.size(); ++i) {
        std::cout << x[i] << ',' << mhat[i] << '\n';
    }
    std::cout << "\n\n";
    for (size_t i = 0; i < y.size(); ++i) {
        std::cout << x[i] << ',' << y[i] << '\n';
    }
    std::cout << "\n\n";
    for (size_t m = 0; m < x.size(); ++m) {
        for (size_t n = 0; n < x.size(); ++n) {
            double scalar_product = 0.0;
            // вычисляем скалярное произведение
            for (size_t i = 0; i < x.size(); ++i) {
                scalar_product += y[i] * psi(i, m, n);
            }
            W[m][n] = scalar_product;
        }
    }
    for (size_t i = 0; i < W.size(); ++i) {
        for (size_t j = 0; j < W[i].size(); ++j) {
            std::cout << i << ',' << j << ',' << W[i][j] << '\n';
        }
    }
    std::cout << "\n\n";
}

int main()
{

    //Test();
    //Task4();
    ////Task5();
    /*FiltrTask1();
    FiltrTask2();
    FiltrTask3();
    FiltrTask4();*/

    ExampleWaveletHAAR();
}
