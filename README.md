
# **Arithmetic Operations on Taylor Series Coefficients**

## Introduction
Higher derivatives have been extracted from Taylor series for the solution of an equation and differential equations. Propagating truncated Taylor series coefficients by evaluation procedure have been used for long time. The d-th derivative can been obtained by recalling the corresponding coefficient [1]. For example, if $f(x)$ and $g(x)$ are both differentiable for n times, $f(x)×g(x)$ is also differentiable for n times. However, using propagating truncated univariate Taylor series would likely result less overflow. Also, the numerical method to find the d-th derivatives can lead to accumulated round off-errors. To avoid mentioned problem arithmetic operations on Taylor series coefficients are suggested, which is a step toward automatic differentiation. 

## Methodology
First, it is necessary to calculate the coefficients of the Taylor series (TS) to generate a function. For any continuously differentiable input path $\text{x(t)}:(-\epsilon,\epsilon)\to \mathbb{R}^n$ with $\epsilon>0$ intermediate variable v by $v(x(t)):(-\epsilon,\epsilon)\to \mathbb{R}^n$ can be defined by the following formula:

$$v(t)=v_0+tv_1+t^2 v_2+\dots+t^d v_d=v(x(t))+o(t^d)$$

Where $v_k$ represents the TS coefficients $(0≤k≤d)$. The coefficients of the function will be calculated in a C++ class and the coefficients will be saved in a vector.

Suppose that $v(t)= φ(u(t),w(t))$ can be obtained either as $v(t)=u(t)+cw(t)$ (with some real constant c) or as $v(t)=u(t)×w(t)$, $v(t)=u(t)/w(t)$, $v(t)=u^2 (t)$, $v(t)=\sqrt{u(t)}$. To obtain the recurrences by applying truncated polynomial arithmetic, we can perform the following operations [1, 2].



$$v_k = u\pm cw = u_k\pm cw_k$$

$$v_k = u × w = \sum_{j=0}^k  u_j w_{k-j}$$

$$v_k = u/w = \frac{u_k-\sum_{j=0}^{k-1}  v_j w_{k-j}}{w_0}$$

$$v_k = u^2 = \sum_{j=0}^k  u_j u_{k-j}$$  

$$v_k = \sqrt{u} =   \frac{u_k-\sum_{j=1}^{k-1}  v_j v_{k-j}}{2v_0}$$  

$$v_k = u^r = \frac{\sum_{j=0}^{k-1} (r(k-j)-j)u_{k-j}v_j}{(ku_0)}$$


The idea is to use the TS coefficients of elementary functions to compute the TS coefficients of a more complex function. For example, the TS coefficients of $f(x) = x$, around $x_0=0$ is equal to {0 , 1}. By using the arithmetic operators mentioned above, it is possible to calculate any polynomial function by this vector. For example, $g(x)=x^2$ is equal to $x×x$.

In case of $u^r$ when $r=0.5$, the first element of vector $u$ cannot be zero (division by zero). For instance, if $f(x)=\sqrt{x}$, the TS representing this function would be an infinite summation of terms.

To test the `ODE_solver()` function, consider the following example: 

$$y' = Q(x) = (f_1(x))-\frac{(f_2(x))}{(f_3(x))}-(f_4(x))×(f_5)^{3.2}(x)+(f_6(x))$$

Where:

$$f_1(x) = 5-x^3+2.5×x+5$$
$$f_2(x) = x^3+3.2-2.5×x^5-1.2+x$$
$$f_3(x) = x^2+1.5$$
$$f_4(x) = 2.5×x^3+x^2+3.5$$
$$f_5(x) = x^2+1.5$$
$$f_6(x) = 4.5×x+8.5$$

The arithmetic operations on $f_1, f_2, f_3, f_4, f_5, f_6$ will result $Q(x)$. With initial value of $y_0 =y(LB=0) = 1$, and the number of divisions of ND=100, The value of $y(UB=0.9)$ will be calculated with an absolute tolerance of AbsTol = $10^{-15}$. The initial guess for all $y_i$ can be considered as equal to zero. All mentioned data besides the number of terms that must be calculated in TS are given in the command line. A new class has been made `TSfinal` to represent $y'$. This class, along with other data, are sent to `ODE_solver()`. The method which has been used to solve the Ode is Euler implicit.

$$y(x_0+(i+1)×h)=y(x_0+ih)+(h/2)(Q(x_0+ih)+Q(x_0+(i+1)×h))$$

Where $h =  (UB-LB)/ND$ is the step size. $Q(x)$ is calculated by `TSfinal` class.

## Code description

There is a single file which contains both the main function and all classes. The constructor `TaylorSeries()` makes the class, while there is a vector named `a` defined in the private. This vector represents the coefficients of the TS.

For any f(x) function, the TS coefficients $(v_k = f^{k}(x_0)/k!)$ are calculated by `tsx` class which represents the TS coefficients of $x$. Since TS is computed around zero, $x_0$ is equal to 0.

Finally, the function `evaluate()` inside the class, computes the value of the given function at any $x=a$. However, the division operator $(/)$  and $\sqrt{}$ is applicable when the $t^k$ converges to zero. Therefor, if $|t|<1$, the formula mentioned above can be used.

To avoid transforming the class to a vector during arithmetic operations on TS coefficients, an operator “[]” is defined, to enhance the time cost.

The input function is given by a txt file named "function.txt", which must be in the same folder. Function `reading_file()` receives the $y'$ as string, and checks the string to make sure it is only a combination of polynomial equations. Functions `PolyN()` and `TOT_poly()` apply the operators by finding the operators in the string.

To make exe file in the folder, following codes should be entered in terminal.

`g++ ped1.cpp -o main`

Other inputs is received from command line, and following structure must be used to send the input to the program.

`./main nTS LB UB y0 ND AbsTol`

Where, nTS is the number of terms needed to be considered in TS. As an example, the following input can be used for the given ODE example and with 120 terms of TS.

`./main 120 0 0.9 1 150 10^-15`

Operators mentioned above $(\pm,/,×,$ and ^$)$ are implemented as functions, in which the output is still a class. At the beginning of each operator a new class `TaylorSeries gg` is made and the size of `a` vector is initialized (the same size of other input). Two `[]` operators are introduced, one of which allows the changes on `a` vector which is necessary for the output `gg[i]`, and the other one is for the inputs (`v` or `w`). To calculate the function, `evaluated()` is used from the class while the inputs are the number of TS terms and the given point $x$ where the function will be calculated at such value. 


## Results
The following table shows the results when ND = 100. A MATLAB code is also developed to compare the efficiency of the results. In MATLAB, the functions representing polynomial functions are used instead of arithmetic operations.

##### Table 1: A comparison of the different number of terms used in TS and MATLAB

|   outputs   |   $nTS = 90$   |   $nTS = 120$   |   MATLAB   |
|-------------|----------------|-----------------|------------|
|   y(0.9)   | -8.22275957406612866 | -8.2227595740661431 | -8.222759574066155 |

It is evident that when $nTS= 90$, it can still cause an accumulated truncation error in ode solver. The reason behind that is calculating TS of a function at a point close to 0.9 which can increase the error due to the operator $/$ and $u^r$ when r is not an integer. This can be problematic especially in ode solver where due to the iterations the error can be accumulated.

It was seen by increasing the nTS to the numbers higher than 120, the accuracy of the decimal numbers was not changed, showing an upper limit for accuracy. However, it is expected to be equal to that of in MATLAB. When the UB was change to 0.5, the value of y(0.5) was equal in both MATLAB. It means the terms with higher power are converging to zero, causing less errors.

## Deficiencies

The first problem is that TS is developed only around $x_0=0.$ Using the operators $/$ and $u^r$ when r in not an integer, can cause noticeable errors. To tackle this issue, it is necessary to develop the general formulation for TS, or calculate the TS for all point between zero to the final point with small intervals.

The function `TOT_poly()` cannot handle is "()" inside another "()". For example, following equation is a valid form.

$$f_1 + (f_2)/(f_3)^r*(f_4)$$

But to read a function similar $z(x)$, based on the number of functions (n), different number of classes are needed.

$$z(x) = (f_1 + (f_2)*((f_3)+((f_4)/(...(f_n)...)))/(f_5))$$

To read this function a `vector<TaylorSeries>` in a loop is necessary, to keep each class in a separate element of the vector.

All $(f_n)^r$ or $x^r$ which are multiplied by a number must be after the number. For example, if $c$ is a constant real number, the valid forms that `TOT_poly()` can read is $c*(f_n)^r$ or $c*x^r$, not $(f_n)^r*c$ or $x^r*c$. `TOT_poly()` cannot handle the error!

For future works, other functions, for example $sin(u), cos(u), log(x), exp(x)$ can be also added as operators.

## Appendix
MATLAB code for ode solver, using MATLAB functions to calculate $Q(x)$ instead of `TaylorSeries()` class defined in C++.

```matlab
clear
clc

aa = 0;
bb = 0.9;
N = 100;
y0 = 1;
abs_Tol = 0.000000000000001;
IG = 0;

f1 = @(x,y) 5-x^3+2.5*x+5-(x^3+3.2-2.5*x^5-1.2+x)/(x^2+1.5)-(2.5*x^3+x^2+3.5)*(x^2+1.5)^3.2+(4.5*x+8.5);

format long

step_size = (bb - aa) / N;
IG1 = IG; % initial guess for all yi
yn1 = IG;
yn0 = y0;

tic
for i=0:N-1

    AT = 100;
    yn1 = IG; % for example when from y0, y1 is calculated, y2 needs to be equal to initial guess
    while (AT > abs_Tol)
        IG1 = yn1; % saving the last iter cal for yn+1
        yn1 = yn0 + step_size * (f1((aa + i * step_size)) + f1((aa + (i + 1) * step_size))) / 2;
        AT = abs(yn1 - IG1); %% until yn1 converges to it's previous iter IG1
    end
    %IG1 = IG;
    yn0 = yn1; % saving the last iter cal for yn+1 in yn, since in next point yn+1 is used as previous point
end

yn_final =yn1;
toc
yn_final

```

## References 
1. Andreas Griewank and Andrea Walther. 2008. Evaluating Derivatives: Principles and Techniques of Algorithmic Differentiation (Second. ed.). Society for Industrial and Applied Mathematics, USA.

2. Ole Stauning. 1996. ENCLOSING SOLUTIONS OF ORDINARY DIFFERENTIAL EQUATIONS : WITH APPLICATIONS. Technical report, Department of Mathematical Modeling, Technical University of Denmark. http://www.imm.dtu.dk/documents/ftp/tr96/tr18_96.abstract.html



