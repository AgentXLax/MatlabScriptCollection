%%
% Nathan Wisla
%
% Numerical Methods (AUMAT 340)
%
% Assignment 3
% 
% March 1, 2016

%% 1. Proving Simpson's Rule
% Simpson's Rule:
%
% $$\int_a^b{f(x)dx} \approx \frac{b-a}{6}(f(x_0 - h) + 4 f(x_0) + f(x_0 + h))=\frac{h}{3}(f(x_0 - h) + 4 f(x_0) + f(x_0 + h))$$
%
%%% a)
% Express the coefficients of $P_2(x)=a_0 + a_1x + a_2x^2$ in terms of $f_0, f_h, f_{-h}$
% 
% $$f(x_0 - h) = a_0 + a_1 (x_0 - h) + a_2 (x_0 - h)^2 $$
% 
% $$f(x_0) = a_0 + a_1 x_0 + a_2 x_0^2 $$
%  
% $$f(x_0 + h) = a_0 + a_1 (x_0 + h) + a_2 (x_0 + h)^2 $$
% 
%%% b) 
% Replace $f(x)$ with $P_2(x)$ in the integral, then integrate the expression.
%%
syms f(x) h c

A = [1 (c - h) (c - h)^2 f(c - h); 1 c (c)^2 f(c) ; 1 (c + h) (c + h)^2 f(c + h)];
A = rref(A);
a0 = A(1,4)
a1 = A(2,4)
a2 = A(3,4)
P2(x) = a0 + a1*x + a2*x^2;
Simp(x) = int(P2,[c-h c+h])

%% 2. Monte-Carlo Method to Estimate pi
% 
% Recall that the area of a circle is $A = \pi r^2$
% Employing Monte-Carlo methods on the unit circle will give us a resulting
% area of $\pi$.
%
% We start by taking a square with side lengths of 2 over the intervals
% $[-1,1]$ for both $x$ and $y$. Then take a set number of random data
% points that lie within that box. If those points land within the
% constraints of the unit circle, we will add them up. If they miss, we
% will ignore them by using the piecewise function
% 
% 
% $$
% C(x,y) = 1, x^2 + y^2 \le 1
% $$
%
% $$
% C(x,y) = 0,  else
% $$
%
% 
% To be more consistent, we could take the average over 1000 Monte-Carlo estimations 
%%
 estPi = zeros(1000,1);
 for i=1:1000
     estPi(i) = monteCarloPi(1000);
 end
 avgPi = mean(estPi)
 absError = abs(pi - avgPi)
 relError = abs((pi - avgPi)/pi)
  
%% 3. Numerical Extrapolation
%
% Given the formula $A = A(h) + Kh^k + \mathcal{O}(h^{k + 1})$,  
% identify $A$ with $f'(x)$ and $A(h)$ with the forward derivative formula
% of $f'(x)$
%
% $$f'(x) = \frac{f(x + h) - f(x)}{h} + Kh^k + \mathcal{O}(h^{k + 1})$$
%
%%% a)
% Find the value of $k$
%
% Here is the equation we were given, rearranged for convenience:
% 
% $$f(x + h) = f(x) + f'(x)h + Kh^k + \mathcal{O}(h^{k + 1})$$
%
% We will compare this equation to a Taylor expansion:
% 
% $$f(x + h) = f(x) + f'(x)h + \frac{f''(x)}{2}h^2 + \mathcal{O}(h^3)$$
% 
% As we can see, the value of $k$ for this extrapolation corresponds to
% $k = 2$.
%
%%% b)
% Now approximate $f'(x)$ using the equation $B(h) = \frac{2^k A(h/2) - A(h)}{2^k - 1}$
%
% Note that $A(h) = \frac{f(x + h) - f(x)}{h}$ and $A = B(h) + \mathcal{O}(h^{k + 1})$
% 
% $$f'(x) = A = B(h) = \frac{2^{2}\frac{f(x + h/2) - f(x)}{h/2} - \frac{f(x + h) - f(x)}{h}}{2^2 - 1}$$
%
% $$f'(x) = \frac{8f(x + h/2) - f(x + h)}{3h}$$
%
% which bears a slight resemblance to the central difference formula.
%
%% 4. Least Squares Fitting
%
% Consider the function $g(x) = x^2 + 1.0$ and its values at 3 distinct
% points: $x_0 = 0,\ x_1 = 1.0,\ x_2 = 2.0$
% 
% Now use least square minimization to fit the function
%
% $$f(x)=ae^{bx}$$ 
% 
% as best as you can to $g(x)$ at those 3 points
%%% a) and b)
% Write down the function $F(a,b)$ to be minimized then define $w = e^b$ and re-write $F(a,b)$ as $G(a,w)$.
% 
% $$F(a,b) = \sum_{i=0}^{2}{(ae^{x_i{\ln{w}}} - g(x_i))^2}$$
%
% $$G(a,w) = \sum_{i=0}^{2}{(aw^{x_i} - g(x_i))^2}$$
%
%%% c)
% write a script to determine the best values for $a$ and $b$ that
% minimize $F$. Also plot $f(x)$ and $g(x)$
%%
g = @(x) x.^2 + 1.0;
xData = [0 1.0 2.0];

F = @(x)sumFunc(x,xData,g(xData)); 
x0 = [0,0];
min = fminsearch(F,x0);
a = min(1)
b = min(2)
f = @(x)(a*exp(b*x));

fplot(f,[-1 2])
hold on
fplot(g,[-1 2])
title('Exponential that best fits 1 + x^2')
legend('y = ae^{bx}','y = 1 + x^2')
xlabel('x')
ylabel('y')

%%
% Is this a linear problem?
%
% $$f(x)=ae^{bx} \Rightarrow \vec{\nabla} f = \vec{0}$$
%
% $$\vec{\nabla} f = \langle{e^{bx}},{ae^{bx}}\rangle \ne \vec{0}$$
% 
% As we can clearly see, the problem cannot be minimized properly because 
% the derivative cannot be set to zero. The problem is not linear.
%% Function Appendix

% Question 2 functions that do Monte-Carlo integration for the unit circle
function f = circleTrue(x,y) % checks to see if the points lie within the unit circle
    if x^2 + y^2 <= 1
       f = 1;
    else 
       f = 0;
    end
end

function f = monteCarloPi(N)
    r = -1 + 2.*rand(N,2); % Generates 2*N random numbers within the bounds [-1,1]
    sum = 0;
    for i = 1:N
       sum = sum + circleTrue(r(i,1),r(i,2));
    end
    
    f = 4/N * sum; % Area of monte carlo box/N times the circleTrue function
end

% Question 4 function
% The vector x is the parameter that will be minimized by fminsearch
function F = sumFunc(x, xData, yData)
    a = x(1);
    b = x(2);
    F = sum((a*exp(b*xData) - yData).^2);
end