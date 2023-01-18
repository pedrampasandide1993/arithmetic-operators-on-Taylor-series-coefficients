#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <string>

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
     * @brief
     * In some cases TS class in made by giving the the size of the vector
     * @param SizeV the size of vector
     */
    TaylorSeries(const uint64_t SizeV)
    {
        make_zero_vec(SizeV);
    }

    /**
     * @brief setting the coefficients of the TS class to the given value
     *
     * @param SizeV Number of TS terms
     * @param value The given value
     */
    void setCoefficient(const uint64_t SizeV, const long double value)
    {
        for (uint64_t i = 0; i < SizeV; i++)
        {
            a[i] = value;
        }
    }
    /**
     * @brief Computing the coefficients of the i-th term for Taylor series
     * @param value The coefficients of a function are given by the input
     */
    void setCoefficient(const vector<long double> &value)
    {
        for (uint64_t i = 0; i < value.size(); i++)
        {
            a.push_back(value[i]);
        }
    }

    /**
     * @brief Using Taylor Series to calculate the value of a function at a given point
     * @param x The point which the value of a function will be calculated
     * @param n The number of derivatives is equal to the number of terms in TS minus 1
     * @param a [i] is the (1/i!)*f^(i)
     * @return long double result the value of the function
     */
    long double evaluate(const long double &x, const uint64_t &n)
    {
        long double result = 0;

        result = a[0] * 1;

        long double x_power_i = 1;

        for (uint64_t i = 1; i <= n; i++)
        {
            x_power_i = x_power_i * x;
            result += a[i] * x_power_i;
        }
        return result;
    }

    /**
     * @brief Initializing the size of vector doing operators on coefficients by operator []
     * Not having the function will result Segmentation error
     * Any v=u*w new TS class (v) will have the same size of u and w, and it is initialized with all zero elements
     * @param index The size of coefficient vector for input(s)
     */
    void make_zero_vec(const uint64_t index)
    {
        a.resize(index);
    }

    /**
     * @brief Finding the size of the coefficient vector
     * Any v=u*w new TS class (v) will have the same size of u and w
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
    vector<long double> a;
};

TaylorSeries operator+(const TaylorSeries &v, const TaylorSeries &w);
TaylorSeries operator+(const TaylorSeries &v, const long double &w);
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
    TaylorSeries gg(v.SIZE_vector());

    for (uint64_t i = 0; i < v.SIZE_vector(); i++)
    {
        gg[i] = v[i] + w[i];
    }

    return gg;
}

TaylorSeries operator+(const TaylorSeries &v, const long double &w)
{
    TaylorSeries gg(v.SIZE_vector());
    gg[0] = v[0] + w;

    for (uint64_t i = 1; i < v.SIZE_vector(); i++)
    {
        gg[i] = v[i];
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
    TaylorSeries gg(v.SIZE_vector());

    for (uint64_t i = 0; i < v.SIZE_vector(); i++)
    {
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
    TaylorSeries gg(v.SIZE_vector());

    for (uint64_t i = 0; i < v.SIZE_vector(); i++)
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
    TaylorSeries gg(w.SIZE_vector());

    for (uint64_t i = 0; i < w.SIZE_vector(); i++)
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

    TaylorSeries gg(v.SIZE_vector());

    gg[0] = v[0] / w[0];
    for (uint64_t i = 1; i < v.SIZE_vector(); i++)
    {
        long double res = 0;
        for (uint64_t j = 1; j <= i; j++)
        {
            res += w[j] * gg[i - j];
        }
        gg[i] = (v[i] - res) / w[0];
    }

    return gg;
}

/**
 * @brief Overloading operator ^ for u = v^w
 * @param v The TS coefficients of the function
 * @param w While const long double would be as the input but only w = 0.5 or 1/2 or 2 would be acceptable when f(x)=x
 * @return TaylorSeries The class representing new function with TS coefficients of u
 */
TaylorSeries operator^(const TaylorSeries &v, const long double &w)
{
    TaylorSeries gg(v.SIZE_vector());

    if (w == 2)
    {
        for (uint64_t i = 0; i < v.SIZE_vector(); i++)
        {
            long double res = 0;
            for (uint64_t j = 0; j <= i; j++)
            {
                res += v[j] * v[i - j];
            }
            gg[i] = res;
        }
    }
    else if (w == 0.5 or w == 1 / 2)
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

            for (uint64_t i = 2; i < v.SIZE_vector(); i++)
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

        TaylorSeries tx0(v.SIZE_vector());
        tx0 = v;
        long double x0 = tx0.evaluate(0, v.SIZE_vector() - 1);
        gg[0] = pow(x0, w);
        if (v.SIZE_vector() >= 2)
        {
            if (gg[0] == 0)
            {
                red();
                cout << "\nError: Division by zero. In f^a(x), the first coefficient of Taylor series for f(x) cannot be zero.\n";
                cout << "\n";
                resetColor();
                exit(-1);
            }

            for (uint64_t i = 1; i < v.SIZE_vector(); i++)
            {
                long double res = 0;
                for (uint64_t j = 0; j <= i - 1; j++)
                {
                    res += (w * (i - j) - j) * v[i - j] * gg[j] / (i * v[0]);
                }
                gg[i] = res;
            }
        }
    }
    return gg;
}

/**
 * @brief Any polynomial equation enclosed with () is sent to this function
 * Any polynomial equation can be represented by the following equation
 * [a_0*x^0] + [a_1*x^1] + [a_2*x^2] ... + [a_d*x^d]
 * each block []_i = a_i*(x^i) is calculated separately from left to the right in PolyN()
 * where a_i can be any real number
 * There is no polynomial equation representing x^a, if 'a' is not an integer. for example x^(0.5)
 * and if a!=2 but it is an integer, x^a is equal to multiplication of x for a times
 * @param LD_b or a_i is the coefficient behind each term, initially it is always +1
 * @param tsx the TS coefficients of f(x)=x
 * @param tsf the TS class to save the TS coefficients of new block
 * @param tsf2 the TS class to save the TS coefficients of summation, tsf2 = tsf2 + tsf
 * @param final_TS_terms the number of TS term needed to be calculated
 * @param str_ped the given string representing a0 + a1*x + x2*x^3 + ... + ad*x^d
 * @return tsf2 the TS class of a0 + a1*x + x2*x^3 + ... + ad*x^d
 */
TaylorSeries PolyN(const TaylorSeries &tsx, const uint64_t &final_TS_terms, const string &str_ped)
{
    TaylorSeries tsf(final_TS_terms);
    TaylorSeries tsf2(final_TS_terms);

    long double LD_b = 1;
    string strNUM;
    for (uint64_t i = 0; i < str_ped.size(); i++)
    {
        string strNUMpow;

        if (str_ped[i] == '+')
        {
            LD_b = 1;
        }
        else if (str_ped[i] == '-')
        {
            LD_b = -1;
        }
        else if (str_ped[i] == '*')
        {
            if (i == 0 || i == str_ped.size() - 1)
            {
                red();
                cout << "\nError: * at the beginning or end of equation.\n";
                cout << "\n";
                resetColor();
                exit(-1);
            }
        }
        else if (str_ped[i] == '^')
        {
            if (i == 0 || i == str_ped.size() - 1)
            {
                red();
                cout << "\nError: ^ at the beginning or end of equation.\n";
                cout << "\n";
                resetColor();
                exit(-1);
            }

            if (isdigit(str_ped[i + 1]) == 0)
            {
                red();
                cout << "\nError: ^ to non-digit value.\n";
                cout << "\n";
                resetColor();
                exit(-1);
            }

            strNUMpow = str_ped[i + 1];
            uint64_t j_pow = 1;
            uint64_t dot_c = 0;

            while (isdigit(str_ped[i + 1 + j_pow]) == 1 || str_ped[i + 1 + j_pow] == '.')
            {
                if (str_ped[i + 1 + j_pow] == '.')
                {
                    dot_c += 1;
                }
                if (dot_c == 2)
                {
                    red();
                    cout << "\nError: ^ to non-digit value. An extra . in the number.\n";
                    cout << "\n";
                    resetColor();
                    exit(-1);
                }
                strNUMpow += str_ped[i + 1 + j_pow];
                j_pow += 1;
            }
            i += (j_pow);
            if (stold(strNUMpow) == stoull(strNUMpow))
            {

                if (stoull(strNUMpow) == 2)
                {
                    tsf = (LD_b * (tsx ^ (stoull(strNUMpow))));
                    LD_b = 1;
                }
                // otherwise x^a = x*x*...x for a times, since in u^r, u0 cannot be zero
                else
                {
                    tsf = tsx;
                    for (uint64_t k = 1; k < stoull(strNUMpow); k++)
                    {
                        tsf = tsf * tsx;
                    }
                    tsf = LD_b * tsf;
                    LD_b = 1;
                }
            }
            else
            {
                red();
                cout << "\nError: x^a cannot be expressed in terms of standard Taylor Series, when a is not a positive integer.\n";
                cout << "The Taylor series of mentioned function is an infinite sum of terms of the form (f^(n)(a))/(n!) * (x-a)^n.\n";
                cout << "\n";
                resetColor();
                exit(-1);
            }

            tsf2 = tsf2 + tsf;
        }
        else if (isdigit(str_ped[i]) == 1)
        {
            strNUM = str_ped[i];
            uint64_t j = 1;
            uint64_t dot_c = 0;
            while (isdigit(str_ped[i + j]) == 1 || str_ped[i + j] == '.')
            {
                if (str_ped[i + j] == '.')
                {
                    dot_c += 1;
                }
                if (dot_c == 2)
                {
                    red();
                    cout << "\nError: ^ to non-digit value. An extra . in the number.\n";
                    cout << "\n";
                    resetColor();
                    exit(-1);
                }

                strNUM += str_ped[i + j];
                j += 1;
            }
            LD_b = LD_b * stold(strNUM);

            // if the one after, i+j+1 is *, the LD_b should be kept, otherwise add it to class, even the end of class
            if (str_ped[i + j] == '-' || str_ped[i + j] == '+' || (i + j - 1) == str_ped.size() - 1)
            {
                tsf2 = tsf2 + LD_b;
                LD_b = 1;
            }

            i += (j - 1);
        }
        else if (str_ped[i] == '.')
        {
            red();
            cout << "\nError: An extra . found in the equation. I must be only between digits.\n";
            cout << "\n";
            resetColor();
            exit(-1);
        }
        else if (str_ped[i] == 'x')
        {
            // if the x is last one or, after that it is + or -, otherwise continue
            if (i == str_ped.size() - 1 || str_ped[i + 1] == '+' || str_ped[i + 1] == '-')
            {
                tsf2 = tsf2 + LD_b * tsx;
                LD_b = 1;
            }
        }

        else
        {
            red();
            cout << "\nError:" << str_ped[i] << "found in the equation.\n";
            cout << "\n";
            resetColor();
            exit(-1);
        }
    }

    return tsf2;
}

/**
 * @brief The inital input string is sent to this function
 * TOT_poly() work like PolyN() but it checks for '(' and ')'
 * And it sends the polynomial equation enclosed with () to PolyN()
 * Based on the next character after ')', including '+','-','*','/','^', it decides to add it to previous TS (ts1)
 * Or keep it in ts2 in case of the one after ')' or the one after digits after '^' is '*' or '/'
 * @param tsx the TS coefficients of f(x)=x
 * @param final_TS_terms the number of TS term needed to be calculated
 * @param str_ped the given string representing the whole
 * @return TSfinal is TS class of a0 + a1*x + x2*x^3 + ... + ad*x^d representing the whole function
 */
TaylorSeries TOT_poly(const TaylorSeries &tsx, const uint64_t &final_TS_terms, const string &str_ped)
{
    TaylorSeries ts1(final_TS_terms);
    ts1.setCoefficient(final_TS_terms, 0);
    TaylorSeries ts2(final_TS_terms);
    ts2.setCoefficient(final_TS_terms, 0);
    TaylorSeries ts3(final_TS_terms);

    TaylorSeries TSfinal(final_TS_terms);

    string strNUMpow;
    long double LD_b = 1;
    string strNUM;
    string strPoly;

    for (uint64_t i = 0; i < str_ped.size(); i++)
    {
        // cout << str_ped[i];
        if (str_ped[i] == '+')
        {
            LD_b = 1;
        }
        else if (str_ped[i] == '-')
        {
            LD_b = -1;
        }

        else if (str_ped[i] == '(')
        {
            uint64_t parentheses = 1;

            uint64_t j = 1;
            while (str_ped[i + j] != ')')
            {
                strPoly += str_ped[i + j];
                if (i + j == str_ped.size())
                {
                    red();
                    cout << "\nError: A missing ).\n";
                    cout << "\n";
                    resetColor();
                    exit(-1);
                }
                j += 1;
            }

            parentheses = parentheses - 1;
            i += j;

            if (i + 1 == str_ped.size())
            {
                // TSfinal = ts1 + ts2; /// ###############################################
                ts1 = PolyN(tsx, final_TS_terms, strPoly);
                ts1 = LD_b * ts1;
                LD_b = 1;
                TSfinal = ts1 + TSfinal;
            }

            else if (str_ped[i + 1] == '+')
            {
                ts1 = PolyN(tsx, final_TS_terms, strPoly);
                ts1 = LD_b * ts1;
                LD_b = 1;
                TSfinal = ts1 + TSfinal;
            }
            else if (str_ped[i + 1] == '-')
            {
                ts1 = PolyN(tsx, final_TS_terms, strPoly);
                ts1 = LD_b * ts1;
                LD_b = 1;
                TSfinal = ts1 + TSfinal;
            }
            else if (str_ped[i + 1] == '*')
            {
                ts2 = PolyN(tsx, final_TS_terms, strPoly);
                ts2 = LD_b * ts2;
                LD_b = 1;
            }
            else if (str_ped[i + 1] == '/')
            {
                ts2 = PolyN(tsx, final_TS_terms, strPoly);
                ts2 = LD_b * ts2;
                LD_b = 1;
            }
            else if (str_ped[i + 1] == '^')
            {
                if (isdigit(str_ped[i + 1 + 1]) == 0)
                {
                    red();
                    cout << "\nError: ^ to non-digit value.\n";
                    cout << "\n";
                    resetColor();
                    exit(-1);
                }

                strNUMpow = str_ped[i + 1 + 1]; // ############# ^5*2
                uint64_t j_pow = 1;
                uint64_t dot_c = 0;

                while (isdigit(str_ped[i + 1 + j_pow + 1]) == 1 || str_ped[i + 1 + j_pow + 1] == '.')
                {
                    // dot counter
                    if (str_ped[i + 1 + j_pow + 1] == '.')
                    {
                        dot_c += 1;
                    }
                    if (dot_c == 2)
                    {
                        red();
                        cout << "\nError: ^ to non-digit value. An extra . in the number.\n";
                        cout << "\n";
                        resetColor();
                        exit(-1);
                    }
                    strNUMpow += str_ped[i + 1 + j_pow + 1];
                    j_pow += 1;
                }
                i += (j_pow + 1); // +1 for first digit
                // x^a is possible only when a is an integer

                if (str_ped[i + 1] == '*' || str_ped[i + 1] == '/')
                {
                    ts2 = PolyN(tsx, final_TS_terms, strPoly);
                    ts2 = (LD_b * (ts2 ^ (stold(strNUMpow))));
                    LD_b = 1; // or will be reassigned in next term by seeing +-
                }
                else if (str_ped[i + 1] == '+' || str_ped[i + 1] == '-')
                {

                    // TSfinal = ts1 + ts2; ///########################################################################
                    ts1 = PolyN(tsx, final_TS_terms, strPoly);
                    ts1 = (LD_b * (ts1 ^ (stold(strNUMpow))));
                    LD_b = 1;
                    TSfinal = ts1 + TSfinal;
                }
                else if (i + 1 == str_ped.size())
                {
                    ts1 = PolyN(tsx, final_TS_terms, strPoly);
                    ts1 = (LD_b * (ts1 ^ (stold(strNUMpow))));
                    LD_b = 1;
                    TSfinal = TSfinal + ts1;
                }

                // x^a.bc check the one after c for +- or / and * to determine for ts2 or ts1
            }

            strPoly = "";
        }
        else if (str_ped[i] == '/' && str_ped[i + 1] == '(')
        {
            if (i == 0 || i == str_ped.size() - 1)
            {
                red();
                cout << "\nError: / at the beginning or end of equation.\n";
                cout << "\n";
                resetColor();
                exit(-1);
            }

            uint64_t parentheses = 1;

            uint64_t j = 1;
            while (str_ped[i + j + 1] != ')')
            {
                strPoly += str_ped[i + j + 1];
                if (i + j + 1 == str_ped.size())
                {
                    red();
                    cout << "\nError: A missing ).\n";
                    cout << "\n";
                    resetColor();
                    exit(-1);
                }
                j += 1;
            }

            parentheses = parentheses - 1;
            i += j + 1;
            ts3 = PolyN(tsx, final_TS_terms, strPoly);

            if (str_ped[i + 1] == '^')
            {
                if (isdigit(str_ped[i + 1 + 1]) == 0)
                {
                    red();
                    cout << "\nError: ^ to non-digit value.\n";
                    cout << "\n";
                    resetColor();
                    exit(-1);
                }

                strNUMpow = str_ped[i + 1 + 1];
                uint64_t j_pow = 1;
                uint64_t dot_c = 0;

                while (isdigit(str_ped[i + 1 + j_pow + 1]) == 1 || str_ped[i + 1 + j_pow + 1] == '.')
                {
                    if (str_ped[i + 1 + j_pow + 1] == '.')
                    {
                        dot_c += 1;
                    }
                    if (dot_c == 2)
                    {
                        red();
                        cout << "\nError: ^ to non-digit value. An extra . in the number.\n";
                        cout << "\n";
                        resetColor();
                        exit(-1);
                    }
                    strNUMpow += str_ped[i + 1 + j_pow + 1];
                    j_pow += 1;
                }
                i += (j_pow + 1);
                ts3 = (LD_b * (ts3 ^ (stold(strNUMpow))));
            }

            ts2 = ts2 / ts3;

            if (i + 1 == str_ped.size())
            {
                TSfinal = TSfinal + ts2;
            }
            else if (str_ped[i + 1] == '-' || str_ped[i + 1] == '+')
            {
                TSfinal = TSfinal + ts2;
            }

            strPoly = "";
        }

        else if (str_ped[i] == '*' && str_ped[i + 1] == '(')
        {
            if (i == 0 || i == str_ped.size() - 1)
            {
                red();
                cout << "\nError: * at the beginning or end of equation.\n";
                cout << "\n";
                resetColor();
                exit(-1);
            }

            uint64_t parentheses = 1;

            uint64_t j = 1;
            while (str_ped[i + j + 1] != ')')
            {
                strPoly += str_ped[i + j + 1];
                if (i + j + 1 == str_ped.size())
                {
                    red();
                    cout << "\nError: A missing ).\n";
                    cout << "\n";
                    resetColor();
                    exit(-1);
                }
                j += 1;
            }

            parentheses = parentheses - 1;
            i += j + 1;
            ts3 = PolyN(tsx, final_TS_terms, strPoly);

            if (str_ped[i + 1] == '^')
            {
                if (isdigit(str_ped[i + 1 + 1]) == 0)
                {
                    red();
                    cout << "\nError: ^ to non-digit value.\n";
                    cout << "\n";
                    resetColor();
                    exit(-1);
                }

                strNUMpow = str_ped[i + 1 + 1];
                uint64_t j_pow = 1;
                uint64_t dot_c = 0;

                while (isdigit(str_ped[i + 1 + j_pow + 1]) == 1 || str_ped[i + 1 + j_pow + 1] == '.')
                {
                    if (str_ped[i + 1 + j_pow + 1] == '.')
                    {
                        dot_c += 1;
                    }
                    if (dot_c == 2)
                    {
                        red();
                        cout << "\nError: ^ to non-digit value. An extra . in the number.\n";
                        cout << "\n";
                        resetColor();
                        exit(-1);
                    }
                    strNUMpow += str_ped[i + 1 + j_pow + 1];
                    j_pow += 1;
                }
                i += (j_pow + 1);
                ts3 = (LD_b * (ts3 ^ (stold(strNUMpow))));
            }

            ts2 = ts2 * ts3;
            strPoly = "";

            if (i + 1 == str_ped.size())
            {
                TSfinal = TSfinal + ts2;
            }
            else if (str_ped[i + 1] == '-' || str_ped[i + 1] == '+')
            {
                TSfinal = TSfinal + ts2;
            }
        }

        // ####################################################################################################

        else if (str_ped[i] == '*')
        {
            if (i == 0 || i == str_ped.size() - 1)
            {
                red();
                cout << "\nError: * at the beginning or end of equation.\n";
                cout << "\n";
                resetColor();
                exit(-1);
            }
            // code
            // keep LD_b
        }
        else if (str_ped[i] == '^')
        {
            if (i == 0 || i == str_ped.size() - 1)
            {
                red();
                cout << "\nError: ^ at the beginning or end of equation.\n";
                cout << "\n";
                resetColor();
                exit(-1);
            }

            //  code
            // checking the number after ^
            if (isdigit(str_ped[i + 1]) == 0)
            {
                red();
                cout << "\nError: ^ to non-digit value.\n";
                cout << "\n";
                resetColor();
                exit(-1);
            }

            strNUMpow = str_ped[i + 1];
            uint64_t j_pow = 1;
            uint64_t dot_c = 0;

            while (isdigit(str_ped[i + 1 + j_pow]) == 1 || str_ped[i + 1 + j_pow] == '.')
            {
                // dot counter
                if (str_ped[i + 1 + j_pow] == '.')
                {
                    dot_c += 1;
                }
                if (dot_c == 2)
                {
                    red();
                    cout << "\nError: ^ to non-digit value. An extra . in the number.\n";
                    cout << "\n";
                    resetColor();
                    exit(-1);
                }
                strNUMpow += str_ped[i + 1 + j_pow];
                j_pow += 1;
            }
            i += (j_pow); // +1 for first digit
            // x^a is possible only when a is an integer
            if (stold(strNUMpow) == stoull(strNUMpow)) // stold() to long double, stoull() to uint64_t
            {
                // using the operator ^2
                if (stoull(strNUMpow) == 2)
                {
                    ts1 = (LD_b * (tsx ^ (stoull(strNUMpow))));
                    LD_b = 1; // or will be reassigned in next term by seeing +-
                    TSfinal = ts1 + TSfinal;
                }
                // otherwise x^a = x*x*...x for a times, since in u^r, u0 cannot be zero
                else
                {
                    ts1 = tsx;
                    for (uint64_t k = 1; k < stoull(strNUMpow); k++)
                    {
                        ts1 = ts1 * tsx;
                    }
                    ts1 = LD_b * ts1;
                    LD_b = 1;
                    TSfinal = ts1 + TSfinal;
                }
            }
            else // base on the definition first element cannot be zero and for f(x)=x the TSc={0,1}
            {
                red();
                cout << "\nError: x^a cannot be expressed in terms of standard Taylor Series, when a is not a positive integer.\n";
                cout << "The Taylor series of mentioned function is an infinite sum of terms of the form (f^(n)(a))/(n!) * (x-a)^n.\n";
                cout << "\n";
                resetColor();
                exit(-1);
            }
        }
        else if (isdigit(str_ped[i]) == 1)
        {
            strNUM = str_ped[i];
            uint64_t j = 1;
            uint64_t dot_c = 0;
            while (isdigit(str_ped[i + j]) == 1 || str_ped[i + j] == '.')
            {
                if (str_ped[i + j] == '.')
                {
                    dot_c += 1;
                }
                if (dot_c == 2)
                {
                    red();
                    cout << "\nError: ^ to non-digit value. An extra . in the number.\n";
                    cout << "\n";
                    resetColor();
                    exit(-1);
                }

                strNUM += str_ped[i + j];
                j += 1;
            }
            LD_b = LD_b * stold(strNUM);

            // if the one after, i+j+1 is *, the LD_b should be kept, otherwise add it to class, even the end of class

            if (str_ped[i + j] == '-' || str_ped[i + j] == '+' || (i + j - 1) == str_ped.size() - 1)
            {
                TSfinal = TSfinal + LD_b;
                LD_b = 1;
            }

            i += (j - 1); // the new i, skip all the str[i+j] saved as a number ########################################i += (j-1)
        }
        else if (str_ped[i] == '.')
        {
            red();
            cout << "\nError: An extra . found in the equation. I must be only between digits.\n";
            cout << "\n";
            resetColor();
            exit(-1);
        }
        else if (str_ped[i] == 'x')
        {
            // if the x is last one or, after that it is + or -, otherwise continue
            if (i == str_ped.size() - 1 || str_ped[i + 1] == '+' || str_ped[i + 1] == '-')
            {
                TSfinal = TSfinal + LD_b * tsx;
                LD_b = 1;
            }
        }

        else
        {
            red();
            cout << "\nError:" << str_ped[i] << "found in the equation.\n";
            cout << "\n";
            resetColor();
            exit(-1);
        }
    }

    return TSfinal;
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
long double ODE_solver(TaylorSeries &ts4, const long double &aa, const long double &bb, const uint64_t &N, const long double &y0, const long double &abs_Tol, const long double &IG, const uint64_t &n_term)
{
    const long double step_size = (bb - aa) / N;
    long double IG1 = IG;
    long double yn1 = IG;
    long double yn0 = y0;

    for (uint64_t i = 0; i < N; i++)
    {
        long double AT = 100;
        yn1 = IG;
        while (AT > abs_Tol)
        {
            IG1 = yn1;
            yn1 = yn0 + step_size * (ts4.evaluate(aa + i * step_size, n_term) + ts4.evaluate(aa + (i + 1) * step_size, n_term)) / 2;

            AT = abs(yn1 - IG1);
        }

        yn0 = yn1;
    }

    return yn1;
}

/**
 * @brief Receiving the filename and reading line by line, first it checks if the string has only valid characters
 * Then it checks for vector pp which is substrings that are not allowed in an equation
 * @param filename
 * @return string
 */
string reading_file(const string &filename)
{
    string line;
    string str_ped;
    ifstream my_file(filename);
    if (my_file.is_open())
    {
        while (getline(my_file, line))
        {
            line.erase(remove(line.begin(), line.end(), ' '), line.end());
            str_ped = str_ped + line;
        }
        my_file.close();
    }
    else
    {
        red();
        cout << "\nError: Unable to open file.\n";
        cout << "\n";
        resetColor();
        exit(-1);
    }

    for (uint64_t i = 0; i < str_ped.size(); i++)
    {
        if (isdigit(str_ped[i]) == 1 || str_ped[i] == '-' || str_ped[i] == '+' || str_ped[i] == '*' || str_ped[i] == '/' || str_ped[i] == '^' || str_ped[i] == ')' || str_ped[i] == '(' || str_ped[i] == 'x' || str_ped[i] == '.')
        {
            continue;
        }
        else
        {
            red();
            cout << "\nError: Not a function. Found '" << str_ped[i] << "' at index " << i << ".\n";
            cout << "\n";
            resetColor();
            exit(-1);
        }
    }

    vector<string> pp = {"()", ")(", "++", "+-", "+*", "+^", "+/", "-+", "--", "-*", "-^", "-/", "*+", "*-", "**", "*^", "*/", "^+", "^-", "^*", "^^", "^/", "/+", "/-", "/*", "/^", "//"};
    string substring;
    for (uint64_t i = 0; i < pp.size(); i++)
    {
        substring = pp[i];

        size_t found = str_ped.find(substring);
        if (found != string::npos)
        {
            red();
            cout << "\nError: Not a function. Found '" << substring << "' at index " << found << ".\n";
            cout << "\n";
            resetColor();
            exit(-1);
        }
    }

    return str_ped;
}

int main(int argc, char *argv[])
{
    ////############################################ INPUTS #################################################
    if (argc > 7)
    {
        red();
        cout << "\nError: Too many input arguments.\n";
        cout << "\n";
        resetColor();
        exit(-1);
    }
    else if (argc < 7)
    {
        red();
        cout << "\nError: Not enough arguments.\n";
        cout << "\n";
        resetColor();
        exit(-1);
    }

    if (stoull(argv[1]) != stold(argv[1]))
    {
        red();
        cout << "\nError: The number of terms in Taylor Series must be an integer.\n";
        cout << "\n";
        resetColor();
        exit(-1);
    }

    if (stoull(argv[5]) != stold(argv[5]))
    {
        red();
        cout << "\nError: The number of divisions in ODE solver must be an integer.\n";
        cout << "\n";
        resetColor();
        exit(-1);
    }

    const uint64_t final_TS_terms = stoull(argv[1]);
    const long double lower_B = stold(argv[2]);
    const long double upper_B = stold(argv[3]);
    const long double y_at_lower_B = stold(argv[4]);
    const uint64_t n_divisions = stoull(argv[5]);
    const long double abs_Tol = stold(argv[6]);

    string str_ped;

    str_ped = reading_file("function.txt");
    cout << "\n";
    cout << "dy/dx = f'(x) = " << str_ped << "\n";

    ////############################################# TS Coeff ###################################################
    const long double inital_guess = 0;

    vector<long double> x;

    x.resize(final_TS_terms, 0);
    x[1] = 1;

    TaylorSeries tsx;
    tsx.setCoefficient(x);
    TaylorSeries TSfinal(final_TS_terms);

    TSfinal = TOT_poly(tsx, final_TS_terms, str_ped);

    ////############################################# ODE SOLVER ###################################################

    timer t;
    t.start();
    long double ODE_res = ODE_solver(TSfinal, lower_B, upper_B, n_divisions, y_at_lower_B, abs_Tol, inital_guess, final_TS_terms - 1);
    t.end();
    cout << "ODE solver took " << t.seconds() << " seconds.\n";

    cout.precision(18);
    cout << "y(" << upper_B << ") = " << ODE_res << "\n\n";

    return 0;
}
