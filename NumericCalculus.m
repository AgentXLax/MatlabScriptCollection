%%
% Nathan Wisla
%
% Numerical Methods (AUMAT 340)
%
% Assignment 2
% 
% February 3, 2016
% 

%% 1. Plotting forward derivative error vs. central difference
% Consider the function $f(x) = e^{x}$ and its derivative $f'(2)$.
%
% a) Write a Matlab code that computes the relative error $\log_{10}{\mathcal{E}}$
%    between the exact value of $f'(2)$ and the forward derivative
%    approximation and the central difference formula, respectively.
%
% b) Plot the results into one graph
%
% c) Can you interpret the gradients?

% a)
errorForward = log10((forwardDerivExp(2) - exp(2)) / exp(2));
errorCentral = log10((centralDifExp(2) - exp(2)) / exp(2));
% b)
x = -12:-1;
plot(x,errorForward,x,errorCentral)
xlabel('log_{10}{h}');
ylabel('log_{10}{E}');
legend('Forward Derivative','Central Difference');

%%
% c)
%
% The slopes to the right of the minimum values on the plot represents the 
% order $\mathcal{O}(h)$ and $\mathcal{O}(h^2)$. Notice how the slope of
% one plot is twice as steep. Since the plotted parameters are logarithms, the
% slope of something $\mathcal{O}(h^2)$ will be 2, and likewise for $\mathcal{O}(h)$
% , the slope will be 1, which is consistent with what was learned in
% class.
%
% The slope to the left of the minimum values represent rounding errors.
%
%% 2. Show the coefficients for the One-sided Derivative at a Boundary Point
% $af(0) + bf(h) + cf(2h) = (a + b + c)f(0) + h(b + 2c)f'(0) + \frac{h^2}{2}(b + 4c)f''(0) + \mathcal{O}(h^3)$
%
% Where we want only $f'(0)$ on the RHS.
%
% Solve for a, b, and c. (System of Equations)
% 
% $$a + b + c = 0$$  
%
% $$h(b + 2c) = 1$$
%
% $$\frac{h^2}{2}(b + 4c) = 0$$
syms h
A = [1 1 1 0; 0 1 2 1/h; 0 1 4 0]
rref(A)

%% 3. Using the forward derivative and the central difference formula
% Use the forward derivative and the central difference formula to compute
% $f'(0)$ for
%
% $f(x) = x + x^2$. Interpret the results.

syms x h
ForwardDeriv = symbolicFD(x)
CentralDif = symbolicCD(x)

%% 
%
% The result for the central difference formula is exact, where the forward
% derivative has an error term of $h$ (The central difference is more accurate)
%% 4. Computing Integrals
% Let
%
% $$f(x) = \frac{1}{x + 1}$$
%
% and consider the integral
%
% $$\int^1_0{f(x)dx}$$
%%
% a) Use the trapezoidal method to estimate the integral
%
% Trapezoidal Method:
% 
% $$\int^b_a{f(x)dx} = (b - a)\frac{f(a) + f(b)}{2}$$

trapMethod = trapezoid(0,1)

%%
% b) Compute the integral exactly and compute the relative error
%
% Use the Fundamental Theorem of Calculus
%
% Relative Error: $\frac{exact - approx.}{exact}$

exactInt = log(2) - log(1) %ln(1) = 0 but it's good to show it
absErrorB = abs(exactInt - trapMethod); %need this for part (e)
relErrorB = abs((exactInt - trapMethod) / exactInt)


%%
% c) Use the composite trapezoidal rule with $n = 2$ to estimate the
% integral. Compute the relative error
%
% Composite Trapezoidal Rule:
%
% $$\int^b_a{f(x)dx} = h[\frac{1}{2}f(x_0) + [\sum_{i=1}^{n-1}{f(x_i)}] + \frac{1}{2}f(x_n)]$$
%
compTrapMethod = compTrap(0,1,2)
absErrorC = abs(exactInt - compTrapMethod);
relErrorC = abs((exactInt - compTrapMethod) / exactInt)
%%
% d) Repeat for $n = 3$

compTrapMethod2 = compTrap(0,1,3)
absErrorD = abs(exactInt - compTrapMethod2);
relErrorD = abs((exactInt - compTrapMethod2) / exactInt)

%%
% e) Use the results for $n = 1, 2, 3$ to estimate how large $n$ has to be
% in order for the relative error to be less than $\frac{1}{1000}$

hvals = [1 1/2 1/3];%values for h from n = 1 to 3
Eh2 = [absErrorB/hvals(1)^2 absErrorC/hvals(2)^2 absErrorD/hvals(3)^2] %these should converge to a certain value that can be used to find n later

% assuming that the third iteration converges closely enough, we can use
% the value Eh2(3) to estimate the value for n.
 

n = ceil(sqrt(Eh2(3) * 10^3))
checkN = log10(abs(log(2) - compTrap(0,1,n)))

%% Function Appendix

function f = forwardDerivExp(x)%displays 12 limit function outputs for the forward derivative.
    listOfVals = zeros(1,12);
    for i = 1:12 %creates steps for a variable h (10^-1 to 10^-12)
        listOfVals(i) = (exp(x + 10^-(13 - i)) - exp(x)) / (10^-(13 - i));
    end    
    f = listOfVals;
end
function f = centralDifExp(x) %displays 12 limit function outputs for the central difference
    listOfVals = zeros(1,12);
    for i = 1:12 %h == 10^-(13 - i)
        listOfVals(i) = (exp(x + 10^-(13 - i)) - exp(x - 10^-(13 - i))) / (2 * 10^-(13 - i));
    end
    f = listOfVals;
end

function f = Q3f(x) %the function used in Q3
    f = x + x^2;
end

function f = symbolicFD(x) %the forward derivative with symbols instead of values
syms h  
limForm = (Q3f(x + h) - Q3f(x)) / h;   
    f = simplify(limForm);
end

function f = symbolicCD(x) %the central derivative with symbols instead of values
    syms h    
    limForm = (Q3f(x + h) - Q3f(x - h)) / (2 * h);   
    f = simplify(limForm);
end

function f = Q4f(x) %the function used in Q4
    f = 1 / (x + 1);
end

function f = trapezoid(a,b) %the trapezoidal method. a and b are endpoints
    f = (b - a) * (Q4f(a) + Q4f(b)) / 2;
end

function f = compTrap(a,b,n) %the composite trapezoidal method. n is the step size
    sum = 0;
    h = (b - a) / n;
    for i = 1:n - 1
        sum = sum + Q4f(a + i * h);
    end
    f = h * (Q4f(a) / 2 + sum + Q4f(b) / 2);    
end

