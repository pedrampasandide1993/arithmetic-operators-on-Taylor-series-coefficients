#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>
using namespace std;

/**
 * @brief Following functions are used to change the colors of errors to red
 * red() Changing the color to red
 * resetColor() Changing the color to default
 */
void red()
{
    printf("\033[1;31m");
}
void resetColor()
{
    printf("\x1b[m");
}

/**
 * @brief Computing the wall clock time
 * The whole timer class is copied from lecture notes from the following ref
 * @ref https://baraksh.com/CSE701/notes.php#creating-and-configuring-your-first-program
 */
class timer
{
public:
    void start()
    {
        start_time = chrono ::steady_clock ::now();
    }
    void end()
    {
        elapsed_time = chrono ::steady_clock ::now() - start_time;
    }
    double seconds() const
    {
        return elapsed_time.count();
    }

private:
    chrono ::time_point<chrono ::steady_clock> start_time = chrono ::steady_clock ::now();
    chrono ::duration<double> elapsed_time = chrono ::duration<double>::zero();
};

class TaylorSeries
{
public:
    /**
     * @brief Constructor for making the class with private variable 'a'
     */
    TaylorSeries() {}

    /**
     * @brief Computing the coefficients of the i-th term for Taylor series
     * @param value The derivatives of function given by user
     */
    void setCoefficient(const vector<long double> &value)
    {
        long double factorial_k = 1;
        for (uint64_t i = 0; i < value.size(); i++)
        {
            a.push_back(value[i] / factorial_k);
            factorial_k *= (long double)(i + 1);
        }
    }

    // Evaluate the series at x
    /**
     * @brief Using Taylor Series to calculate the value of a function at a given point
     * @param x The point which the value of a function will be calculated
     * @param n The number of derivatives is equal to the number of terms in TS minus one
     * @return long double result The value of the function
     */
    long double evaluate(const long double &x, const uint64_t &n)
    {
        long double result = 0;

        // first term i=0 which 0!=1, with this algorithm using pow() function and factorial is not necessary
        // x^0 = 1 (first 1), 0!=1 (second i),
        result = a[0] * 1;

        long double x_power_i = 1;

        for (uint64_t i = 1; i <= n; i++)
        {
            x_power_i = x_power_i * x;
            result += a[i] * x_power_i; // a[i] is the (1/i!)*f^(i) i-derivative of function f
        }
        return result;
    }

    /**
     * @brief Initializing the size of vector doing operators on coefficients by operator []
     * Not having the function will result Segmentation error
     * @param index The size of coefficient vector for input(s)
     */
    void make_zero_vec(const uint64_t index)
    {
        a.resize(index);
    }

    /**
     * @brief Finding the size of a vector, which is used in previous function
     * @return uint64_t
     */
    uint64_t SIZE_vector() const
    {
        return a.size();
    }

    /**
     * @brief Overloaded operator [] to access vector elements
     * Does not allow the modification on elements
     * @param index The index of element
     * @return const long double& The value of the element
     */
    const long double &operator[](uint64_t index) const
    {
        return a[index];
    }

    /**
     * @brief Overloaded operator [] to access vector elements
     * Allows the modification on elements
     * @param index The index of element
     * @return long double& The new value of the element
     */
    long double &operator[](uint64_t index)
    {
        return a[index];
    }

private:
    /**
     * @brief The TS coefficients of a function
     */
    vector<long double> a; //(1/i!)*f^(i) is i-order derivative of function f
};

/**
 * @brief A function given by the user
 * @param x The point where the function will be calculated
 * @return long double
 */
long double func1(const long double x)
{
    return exp(-x) * x + 3;
}

/**
 * @brief A function given by the user
 * @param x The point where the function will be calculated
 * @return long double
 */
long double func2(const long double x)
{
    return exp(x) + x * x + 1.5;
}

TaylorSeries operator+(const TaylorSeries &v, const TaylorSeries &w);
TaylorSeries operator-(const TaylorSeries &v, const TaylorSeries &w);
TaylorSeries operator*(const TaylorSeries &v, const TaylorSeries &w);
TaylorSeries operator*(const long double &v, const TaylorSeries &w);
TaylorSeries operator*(const TaylorSeries &w, const long double &v);
TaylorSeries operator/(const TaylorSeries &v, const TaylorSeries &w);
TaylorSeries operator^(const TaylorSeries &v, const long double &w);

/**
 * @brief Overloading operator + for u = v+w
 * @param v The TS coefficients of first function
 * @param w The TS coefficients of second function
 * @return TaylorSeries The class representing new function with TS coefficients of u
 */
TaylorSeries operator+(const TaylorSeries &v, const TaylorSeries &w)
{
    TaylorSeries gg;
    gg.make_zero_vec(v.SIZE_vector());

    for (uint64_t i = 0; i < v.SIZE_vector(); i++) // size should be the size of a vector
    {
        // cout << v[i] + w[i];
        gg[i] = v[i] + w[i];
    }

    return gg;
}

/**
 * @brief Overloading operator - for u = v-w
 * @param v The TS coefficients of first function
 * @param w The TS coefficients of second function
 * @return TaylorSeries The class representing new function with TS coefficients of u
 */
TaylorSeries operator-(const TaylorSeries &v, const TaylorSeries &w)
{
    TaylorSeries gg;
    gg.make_zero_vec(v.SIZE_vector());

    for (uint64_t i = 0; i < v.SIZE_vector(); i++) // size should be the size of a vector
    {
        // cout << v[i] + w[i];
        gg[i] = v[i] - w[i];
    }

    return gg;
}

/**
 * @brief Overloading operator * for u = v*w
 * @param v The TS coefficients of first function
 * @param w The TS coefficients of second function
 * @return TaylorSeries The class representing new function with TS coefficients of u
 */
TaylorSeries operator*(const TaylorSeries &v, const TaylorSeries &w)
{
    TaylorSeries gg;
    gg.make_zero_vec(v.SIZE_vector());

    for (uint64_t i = 0; i < v.SIZE_vector(); i++) // size should be the size of a vector
    {
        long double res = 0;
        for (uint64_t j = 0; j <= i; j++)
        {
            res += v[j] * w[i - j];
        }
        gg[i] = res;
    }

    return gg;
}

/**
 * @brief Overloading operator * for u = v*w
 * @param v The scalar number
 * @param w The TS coefficients of the function
 * @return TaylorSeries The class representing new function with TS coefficients of u
 */
TaylorSeries operator*(const long double &v, const TaylorSeries &w)
{
    TaylorSeries gg;
    gg.make_zero_vec(w.SIZE_vector());

    for (uint64_t i = 0; i < w.SIZE_vector(); i++) // size should be the size of a vector
    {
        gg[i] = v * w[i];
    }

    return gg;
}

/**
 * @brief Overloading operator * for u = v*w
 * @param v The scalar number
 * @param w The TS coefficients of the function
 * @return TaylorSeries The class representing new function with TS coefficients of u
 */
TaylorSeries operator*(const TaylorSeries &w, const long double &v)
{
    return v * w;
}

/**
 * @brief Overloading operator / for u = v/w
 * @param v The TS coefficients of first function
 * @param w The TS coefficients of second function
 * @return TaylorSeries The class representing new function with TS coefficients of u
 */
TaylorSeries operator/(const TaylorSeries &v, const TaylorSeries &w)
{
    if (w[0] == 0)
    {
        red();
        cout << "\nError: Division by zero. In f(x)/g(x), the first coefficient of Taylor series for g(x) cannot be zero.\n";
        cout << "\n";
        resetColor();
        exit(-1);
    }

    TaylorSeries gg;
    gg.make_zero_vec(v.SIZE_vector());

    gg[0] = v[0] / w[0];
    for (uint64_t i = 1; i < v.SIZE_vector(); i++) // size should be the size of a vector
    {
        long double res = 0;
        for (uint64_t j = 1; j <= i; j++)
        {
            res += w[j] * gg[i - j];
        }
        gg[i] = (v[i] - res) / w[0];
    }

    return gg; // check the size of gg.a
}

/**
 * @brief Overloading operator ^ for u = v^w
 * @param v The TS coefficients of the function
 * @param w While const long double would be in inputs but only w = 0.5 or 1/2 or 2 would be acceptable
 * @return TaylorSeries The class representing new function with TS coefficients of u
 */
TaylorSeries operator^(const TaylorSeries &v, const long double &w)
{
    TaylorSeries gg;
    gg.make_zero_vec(v.SIZE_vector());

    if (w == 2)
    {
        for (uint64_t i = 0; i < v.SIZE_vector(); i++) // size should be the size of a vector
        {
            long double res = 0;
            for (uint64_t j = 0; j <= i; j++)
            {
                res += v[j] * v[i - j];
            }
            gg[i] = res;
        }
    }
    else if (w == 0.5 or w == 1 / 2) // 1/2?
    {
        gg[0] = sqrt(v[0]);
        if (v.SIZE_vector() >= 2)
        {
            if (gg[0] == 0)
            {
                red();
                cout << "\nError: Division by zero. In f^0.5(x), the first coefficient of Taylor series for f(x) cannot be zero.\n";
                cout << "\n";
                resetColor();
                exit(-1);
            }

            gg[1] = v[1] / (2 * gg[0]);

            for (uint64_t i = 2; i < v.SIZE_vector(); i++) // size should be the size of a vector
            {
                long double res = 0;
                for (uint64_t j = 1; j <= i - 1; j++)
                {
                    res += gg[j] * gg[i - j];
                }
                gg[i] = (v[i] - res) / (2 * gg[0]);
            }
        }
    }
    else
    {
        red();
        cout << "\nError: The power can be either 2 or 0.5.\n";
        cout << "\n";
        resetColor();
        exit(-1);
    }
    return gg;
}

/**
 * @brief Solving the ODE by implicit Euler's Method, y(bb)=?
 * @param ts4 The class with TS coefficients representing a function given by the problem
 * @param aa The lower bound of span given by the problem
 * @param bb The upper bound of span given by the problem
 * @param N The number of divisions given by the problem
 * @param y0 The initial value y(aa)=y0 given by the problem
 * @param abs_Tol The Absolute Tolerance allowed between two iterations given by the problem
 * @param IG The initial guess for all points given by the problem
 * @param n_term The number of terms used in TS class
 * (n_term is the vector size of derivative = number of terms -1)
 * @return long double The value of y(bb)
 */
long double ODE_solver(TaylorSeries ts4, const long double aa, const long double bb, const uint64_t N, const long double y0, const long double abs_Tol, const long double IG, const uint64_t n_term)
{
    const long double stpsize = (bb - aa) / N;
    long double IG1 = IG; // inital guess for all yi
    long double yn1 = IG;
    long double yn0 = y0;

    for (uint64_t i = 0; i < N; i++)
    {
        long double AT = 100;
        yn1 = IG; // for example when from y0, y1 is calculated, y2 needs to be equal to initial guess
        while (AT > abs_Tol)
        {
            IG1 = yn1; // saving the last iter cal for yn+1
            yn1 = yn0 + stpsize * (ts4.evaluate(aa + i * stpsize, n_term) + ts4.evaluate(aa + (i + 1) * stpsize, n_term)) / 2;

            AT = abs(yn1 - IG1); /// until yn1 converges to it's previous iter IG1
        }
        // IG1 = IG;
        yn0 = yn1; // saving the last iter cal for yn+1 in yn, since in next point yn+1 is used as previous point
    }

    return yn1;
}

int main()
{
    ////############################################# INPUTS #################################################

    /**
     * @brief All the following inputs are given by the user
     *
     * const long double aa = 0.5     Computing TS at x=aa, or ODE
     * vector<long double> f_diff1p   The derivatives of f1 at x=0
     * vector<long double> f_diff2p   The derivatives of f2 at x=0
     *
     * const long double lower_B      The lower bound used in ODE
     * const long double upper_B      The upper bound used in ODE
     * const uint64_t n_divisions     The number of divisions in span [lower_B upper_B]
     * const long double y_at_lower_B In IVP ODE y(x=lower_B)=y_at_lower_B
     * const long double abs_Tol      The absolute tolerance error allowed in implicit Euler's method
     * const long double inital_guess The initial guess for all yi in the span
     */
    const long double aa = 0.5;
    const long double lower_B = 1;
    const long double upper_B = 5;
    const uint64_t n_divisions = 100;
    const long double y_at_lower_B = 7.1651;
    const long double abs_Tol = 0.000000000000001;
    const long double inital_guess = 0;

    /**
     * @brief The alternative derivative for function func1 =  x*exp(-x) + 3 at x=0 with different number of terms
     */
    // vector<long double> f_diff1p = {3, 1, -2, 3, -4, 5, -6, 7, -8, 9, -10}; // these are f^0(x=0) to f^10(x=0)
    // vector<long double> f_diff1p = {3, 1, -2, 3, -4, 5, -6, 7, -8, 9, -10, 11, -12, 13, -14, 15, -16, 17, -18, 19, -20, 21, -22, 23}; // these are f^0(x=0) to f^10(x=0)
    vector<long double> f_diff1p = {3, 1, -2, 3, -4, 5, -6, 7, -8, 9, -10, 11, -12, 13, -14, 15, -16, 17, -18, 19, -20, 21, -22, 23, -24, 25, -26, 27, -28, 29, -30, 31, -32, 33, -34, 35}; // these are f^0(x=0) to f^10(x=0)

    /**
     * @brief The alternative derivative for function func2 =  exp(x) + x * x + 1.5 at x=0 with different number of terms
     */
    // vector<long double> f_diff2p = {2.5, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1}; // these are f^0(x=0) to f^10(x=0)
    // vector<long double> f_diff2p = {2.5, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // these are f^0(x=0) to f^10(x=0)
    vector<long double> f_diff2p = {2.5, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // these are f^0(x=0) to f^10(x=0)

    ////######################################################################################################
    const uint64_t n_term = f_diff1p.size() - 1;
    const uint64_t n_term2 = f_diff2p.size() - 1;
    /**
     * @brief n_term is the number of terms except in TS except the first term f(x0), The number of terms in both functions must be equal
     */
    if (n_term2 != n_term)
    {
        red();
        cout << "Error: The number of terms considered in both functions must be equal.\n";
        cout << "\n";
        resetColor();
        exit(-1);
    }

    TaylorSeries ts1;
    ts1.setCoefficient(f_diff1p);

    TaylorSeries ts2;
    ts2.setCoefficient(f_diff2p);

    cout.precision(18);

    cout << "using TS class f(x=a) = " << ts1.evaluate(aa, n_term) << "\n"; // Outputs 2.5 (1 + 1 + 1/2 + 1/6)
    cout << "Absolute Error = " << abs(func1(aa) - ts1.evaluate(aa, n_term)) << "\n\n";

    cout << "using TS class g(x=a) = " << ts2.evaluate(aa, n_term) << endl; // Outputs 2.5 (1 + 1 + 1/2 + 1/6)
    cout << "Absolute Error = " << abs(func2(aa) - ts2.evaluate(aa, n_term)) << "\n\n";

    cout << "using TS class f(x=a) + g(x=a) = " << (ts1 + ts2).evaluate(aa, n_term) << "\n";
    cout << "Absolute Error = " << abs(func1(aa) + func2(aa) - (ts1 + ts2).evaluate(aa, n_term)) << "\n\n";

    cout << "using TS class f(x=a) - g(x=a) = " << (ts1 - ts2).evaluate(aa, n_term) << "\n";
    cout << "Absolute Error = " << abs(func1(aa) - func2(aa) - (ts1 - ts2).evaluate(aa, n_term)) << "\n\n";

    cout << "using TS class f(x=a) * g(x=a) = " << (ts1 * ts2).evaluate(aa, n_term) << "\n";
    cout << "Absolute Error = " << abs(func1(aa) * func2(aa) - (ts1 * ts2).evaluate(aa, n_term)) << "\n\n";

    cout << "using TS class 2 * g(x=a) = " << (2 * ts2).evaluate(aa, n_term) << "\n";
    cout << "Absolute Error = " << abs(2 * func2(aa) - (2 * ts2).evaluate(aa, n_term)) << "\n\n";

    cout << "using TS class g(x=a) *2 = " << (ts2 * 2).evaluate(aa, n_term) << "\n";
    cout << "Absolute Error = " << abs(func2(aa) * 2 - (ts2 * 2).evaluate(aa, n_term)) << "\n\n";

    cout << "using TS class f(x=a) / g(x=a) = " << (ts1 / ts2).evaluate(aa, n_term) << "\n";
    cout << "Absolute Error = " << abs((func1(aa) / func2(aa)) - (ts1 / ts2).evaluate(aa, n_term)) << "\n\n";

    cout << "using TS class f^2(x=a) = " << (ts1 ^ 2).evaluate(aa, n_term) << "\n";
    cout << "Absolute Error = " << abs(pow(func1(aa), 2) - (ts1 ^ 2).evaluate(aa, n_term)) << "\n\n";

    cout << "using TS class f^0.5(x=a) = " << (ts1 ^ 0.5).evaluate(aa, n_term) << "\n";
    cout << "Absolute Error = " << abs(sqrt(func1(aa)) - (ts1 ^ 0.5).evaluate(aa, n_term)) << "\n\n";

    cout << "using TS class f(x=a) + 2 * g(x=a) = " << (ts1 + 2 * ts2).evaluate(aa, n_term) << "\n";
    cout << "Absolute Error = " << abs(func1(aa) + 2 * func2(aa) - (ts1 + 2 * ts2).evaluate(aa, n_term)) << "\n\n";

    TaylorSeries ts4;
    ts4 = (ts1 * ts2);

    timer t;
    t.start();
    long double ODE_res = ODE_solver(ts4, lower_B, upper_B, n_divisions, y_at_lower_B, abs_Tol, inital_guess, n_term);
    t.end();
    cout << "ODE solver took " << t.seconds() << " seconds.\n";

    cout << "y(" << upper_B << ") = " << ODE_res << "\n";
}
