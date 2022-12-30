
# **Arithmetic Operations on Taylor Series Coefficients for Different Functions**

## Introduction
First, we need to calculate the coefficients of the Taylor series (TS) to generate a function. For any continuously differentiable input path $\text{x(t)}:(-\epsilon,\epsilon)\to \mathbb{R}^n$ with $\epsilon>0$ intermediate variable v by $v(x(t)):(-\epsilon,\epsilon)\to \mathbb{R}^n$ can be defined by the following formula:
$$v(t)=v_0+tv_1+t^2 v_2+\dots+t^d v_d=v(x(t))+o(t^d)$$
Where $v_k$ represents the TS coefficients $(0≤k≤d)$. The coefficients of the function will be calculated in a C++ class and the coefficients will be saved in a vector.

Suppose that $v(t)= φ(u(t),w(t))$ can be obtained either as $v(t)=u(t)+cw(t)$ (with some real constant c) or as $v(t)=u(t)×w(t)$, $v(t)=u(t)/w(t)$, $v(t)=u^2 (t)$, $v(t)=\sqrt{u(t)}$. To obtain the recurrences by applying truncated polynomial arithmetic, we can perform the following operations.


##### Table 1: Taylor Series Propagation through Arithmetic Operations
| $v=$       | recurrence for $k=1…d$ |
|----------  |----------|
| $u\pm cw$  | $v_k = u_k\pm cw_k$   |
| $u × w$    | $v_k=\sum_{j=0}^k  u_j w_{k-j}$  |
| $u/w$      | $v_k = [u_k-\sum_{j=0}^{k-1}  v_j w_{k-j}]/w_0$         |
| $u^2$      | $v_k=\sum_{j=0}^k  u_j u_{k-j}$   |
| $\sqrt{u}$ | $v_k = [u_k-\sum_{j=1}^{k-1}  v_j v_{k-j}]/2v_0$   |

The purpose of this project is to reduce the runtime by conducting arithmetic operations on coefficients before calculating the real value of a given function.

## Code description

There is a single file which contains both the main function and class. The constructor `TaylorSeries()` makes the class, while there is a vector named `a` defined in the private. This vector represents the coefficients of the TS.

For any f(x) function, the coefficients of TS ($v_k = f^{k}(x_0)/k!$) are calculated by `void setCoefficient()` inside the class. The input of setCoefficient(), which is given by the user, is a vector of derivatives of the function $(f^{(k)} (x_0))$ where $k=0…d$. Since TS is computed around zero, $x_0$ is equal to 0.

Finally, the function `evaluate()` inside the class, computes the value of the given function at any $x=a$.

To avoid transforming the class to a vector during arithmetic operations on TS coefficients, an operator “[]” is defined, to enhance the time cost.

Functions $f_1(x)$ and $f_2(x)$ are given function by the user as `func1()` and `func2()`, respectively. The program intends to implement arithmetic operations on their coefficients. As an example the following functions have been chosen.

$$f_1(x)=x\times \exp(-x)+3$$
$$f_2(x)=\exp(x)+ x^2+1$$

`f_diff1p` and `f_diff2p` are derivatives of $f_1(x)$ and $f_2(x)$ at $x=0$ given by the user inside the main function (Section INPUTS). As mentioned above, the `setCoefficient()` function takes these coefficients and calculates the TS coefficients.

Operators mentioned in Table 1 $(\pm,/,×,$ and ^$)$ are implemented as functions, in which the output is still the class. At the beginning of each operator a new class `TaylorSeries gg` is made and the size of `a` vector is initialized (the same size of other input). Two `[]` operators are introduced, one of which allows the changes on `a` vector which is necessary for the output `gg[i]`, and the other one is for the inputs (`v` or `w`) which are the TS coefficients of $f_1$ and $f_2$. To calculate the function, `evaluated()` is used from the class while the inputs are the number of TS terms and the given point x where the function will be calculated at such value. 




## Results
The operators were performed for functions func1 and func2 and results are shown in Table 2.

##### Table 2: The result of arithmetic operations on $f_1=f$ and $f_2=g$
| `TaylorSeries()` at $a=0.5$ | Absolute Tol $d=10$ | Absolute Tol $d=10$ |
|----------     |----------              |----------|
| $f(a)$        | 1.28686365816894e-10   | 0   |
| $g(a)$        | 1.27624887654582e-11   | 0  |
| $f(a)+g(a)$   | 1.41448854365511e-10   | 4.33680868994202e-19|
| $f(a)-g(a)$   | 1.1592387675328e-10    | 1.49077798716757e-19|
| $f(a)×g(a)$   | 1.1700852542082e-08    | 8.67361737988404e-19|
| $2×g(a)$      | 2.55249775309163e-11   | 0   |
| $g(a)×2$      | 2.55249775309163e-11   | 0   |
| $f(a)/g(a)$   | 1.20853807229532e-05   | 1.50162000889242e-17|
| $f^2 (a)$     | 6.25046132337291e-07   | 1.73472347597681e-18|
| $f^{0.5} (a)$ | 4.28735177033699e-06   | 5.76166718502247e-15|
| $f(a)+2×g(a)$ | 1.54211342480448e-10   | 8.67361737988404e-19|

The Table 2 shows the absolute tolerance $(|real-TS class|)$ when the number of terms considered in TS is 11 and 36, in which the error might be negligible. In the case of $f^{0.5} (a)$ and $f(a)/g(a)$ the error is higher than the other operators. The reason for this is that the the formulation in Table 1 for these two operators are applicable when the $t^k$ converges to zero. Therefor, if $|t|<1$, the formula mentioned above can be used.

To test the `ODE_solver()` function, consider the following example: 

$$y' = Q(x) = f_1 × f_2 = \left( x × \exp(-x) + 3 \right) \left( \exp(x) + x^2 + 1.5 \right)$$

With initial value of $y(a=1) = 7.1651$, and the number of divisions of N=100. The value of $y(b=5)$ will be calculated with an Absolute Tolerance of $10^{-15}$. The initial guess for all yi can be considered as equal to zero. All mentioned data is given in INPUTS section of the main function. A new class has been made `ts4 = ts1*ts2` to represent $y'$. This class, along with other data, are sent to `ODE_solver()`. The method which has been used to solve the Ode is Euler implicit.

$$y(x_0+(i+1)×h)=y(x_0+ih)+(h/2)(Q(x_0+ih)+Q(x_0+(i+1)×h))$$

While $h =  (b-a)/N$   is the step size. $Q(x)$ is calculated by `ts4` class.

##### Table 3: A comparison of the different number of terms used in TS and MATLAB

| outputs | when $d = 10$ | when $d = 35$| MATLAB |
|----------|----------|----------|----------|
| y(5)   | 265.577195222504461   | 603.649972173915477   | 603.649972173915277   |
| Wall Time (sec)   | $≈0.0002$   | $≈0.0002$   | $≈0.0009$   |

It is evident that when $d= 10$, it can cause a significant accumulated truncation error in ode solver. The MATLAB code used here is mentioned in Appendix. Notably, using a higher precision in MATLAB by vpa() function the time cost increases up to ≈ 2.5 seconds.

## Deficiencies
The Taylor series function is defined by the following formulation. 

$$f(x) = f(x_0) + (x - x_0) f'(x_0) / 1! + (x - x_0)^2 f''(x_0) / 2! + ⋯ + (x - x_0)^d f^d(x_0) / d!$$

The coefficients of the TS $f^d (x_0)/d!$ need to be calculated by differentiations. If we use numerical methods to do this, for example forward differentiation:

$$f^d(x_0)=(\sum_{i=0}^n  (-1)^{n-i} \binom{n}{k} f(x+ih))/h^d$$

After the first derivative the result would not be accurate due to the accumulation of round off-errors.

The second method is automatic differentiation in C++ , like FADBAD++. However, in this project the derivatives are given by the user.

The third solution is using the coefficients of elementary functions to compute the TS coefficients of a more complex function. For example, the TS coefficients of $f(x) = x$, around $x0=0$ is equal to (0 , 1). By using the arithmetic operators mentioned in Table 1 we can calculate any polynomial function by this vector. For example, $g(x)=x^2$ is equal to $x×x$.

The second problem is that TS is developed only around $x0=0.$ Using the operators $/$ and $\sqrt{}$ can cause noticeable errors. To tackle this issue, it is necessary to develop the general formulation for TS, or calculate the TS for all point between zero to the final point with small intervals.

## Appendix
MATLAB code for ode solver, using MATLAB functions to calculate $Q(x)$ instead of `TaylorSeries()` class defined in C++.


```matlab
clear
clc

aa = 1;
bb = 5;
N = 100;
y0 = 7.1651;
abs_Tol = 0.000000000000001;
IG = 0;

f1 = @(x) (x * exp(-x) +3)*(exp(x) + x * x + 1.5);

format long

stpsize = (bb - aa) / N;
IG1 = IG; % initial guess for all yi
yn1 = IG;
yn0 = y0;

tStart = vpa(cputime);
tic
for i=0:N-1

    AT = 100;
    yn1 = IG; % for example when from y0, y1 is calculated, y2 needs to be equal to initial guess
    while (AT > abs_Tol)
        IG1 = yn1; % saving the last iter cal for yn+1
        yn1 = yn0 + stpsize * (f1((aa + i * stpsize)) + f1((aa + (i + 1) * stpsize))) / 2;
        AT = abs(yn1 - IG1); %% until yn1 converges to it's previous iter IG1
    end
    %IG1 = IG;
    yn0 = yn1; % saving the last iter cal for yn+1 in yn, since in next point yn+1 is used as previous point
end

ynfinal =yn1;
toc
ynfinal
tEnd = vpa(cputime) - vpa(tStart)


```




