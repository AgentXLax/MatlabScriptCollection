%% Assignment 5
% Nathan Wisla
%
% Numerical Methods (AUMAT 340)
% 
% April 5, 2016

%% 1) Newton's Method for a system of equations
% Find the roots of the system 
% 
% $$F(x_1,x_2) = \langle f_1(x_1,x_2)\rangle = \langle e^{2x_1-x_2}-x_1,x_1^2-x_2 \rangle$$
%
% By solving the parabolic equation $x_1^2-x_2=0$, and substituting the result
% into the exponential, we have reduced the parameters to make a guess.
% The exponential has a local maximum around $x=1$, with no roots when 
% $x_1<1$ so we must choose an $x_1>1$. For the parabolic function, we can 
% make any guess for $x_2$ because convergence will be stable on a slanted 
% parabola.
%%
syms x y a b
v = [x,y];
f = @(x) [exp(2*x(1)-x(2))-x(1); x(1)^2-x(2)];
J = @(x) jacobian(f(x));
Jsub = @(b) eval(subs(subs(J(v),b(1)),b(2)));

x0 = [1.5;0];

newtonRecursive = newton2DRW(f,x0,10^-6)
newtonIterative = newton2DI(f,x0,10^-6)

%%
% As we can see (bottom of document again...), Newton's method converges.



%% 2) Newton's Method for Extrema
% Consider the function $F(x_1,x_2)=-\ln{1-x_1-x_2}-\ln{x_1}-\ln{x_2}$
% Find the minimum of this function
% First, we will look at the contour plot to visualize where the minimum of
% the function might be. Using the surface plot made it look pretty, so I
% used that instead of the contour method.
syms x y a b
v = [x,y];
F = @(x) -log(1 - x(1) -x(2)) - log(x(1)) - log(x(2));

fsurf(F(v),[0.2 0.4])
view(2); %viewing angle is the x,y plane
%% 
% As we can see, it appears that the function has a minimum around the
% point $(0.33,0.33)$.
x0 = [0.85;0.05];

newtonMaxMin = newton2DMaxMin(F,x0,10^-6)

%%
% The solution turns out to be exactly at that point.
%% 3) 4th Order Runge-Kutta Method
% Solve the differential equation $y'(t) = y^2 \sin{t}, t\in [0,2],  y(0)=0.1$
%
% Runge-Kutta: $y(i+1) = y(i) + 1/6(K1+2*K2+2*K3+K4)$

y0 = 0.1;
% approximate the differential equation
f = @(t,y) y^2 + sin(t);

h=0.1;
dy = rungeKutta4(f,y0,h);
[small big] = size(dy);
while abs(dy(big) - dy(big - 1)) > 10^-6
    disp('checking step size:');
    disp(h);
    dy = rungeKutta4(f,y0,h);
    h = h/10;
    [small big] = size(dy);
end
%%
% Due to time complexity, it takes too long to get proper x-scaling on the
% graph. However, the range is correct with improper orders of magnitude.
disp('step size for 10^-6:');
disp(h);
plot(dy);
%% 4) Implicit Euler Method
% If we want to solve $y'(t)=y\cos{y}, t\in[0,1], y(0)=1$ with the implicit
% Euler method, we need to iterate with the following formula
%
% $$y_{i+1}=y_{i} + hy'_{i+1}$$
%
% $$y_{i+1}=y_{i} + hy_{i+1}\cos{y_{i+1}}$$
%
% $$y_{1}= 1 + hy_{1}\cos{y_{1}}$$
% 
% as we can see, this equation is not linear, so if we want to solve this
% equation, we would have to go about it by using Newton's method on the non-
% linear derivative. A good guess for the root would be $x\in[-1,1]$

%% 5) Boundary-Value Problems
% Solve the BVP with a step size $h=0.1$
% 
% $$u''(x)+u'(x)+u(x)=0,u(0)=1,u(1)=2$$
h = 0.1;
xl = 0;
xr = 1;
u1 = 1;
uN = 2;
N = (xr - xl)/h;
A = zeros(N-1,N-1);
b = zeros(N-1,1);
x = zeros(N+1,1);
u = zeros(N+1,1);

b(1) = 1*(h-2);
b(N-1) = -2*(h+2);
A(1,1) = 2*h^2 - 4;
A(1,2) = (2 + h);
for i = 2:(N-2)
    A(i,i-1) = 2 - h;
    A(i,i) = 2*h^2 - 4;
    A(i,i+1) = 2 + h;
end
A(N-1,N-2) = (2 - h);
A(N-1,N-1) = 2*h^2 - 4;

v = A\b;

x(1) = xl;
x(N+1) = xr;
u(1) = u1;
u(N+1) = uN;
for j = 2:N
    x(j)=(j-1)*h;
    u(j)=v(j-1);
end
plot(x,u)
%%
% BONUS: How would the problem change if we replace $u'(x)$ with
% $u(x)u'(x)$?
%
% The iterative equation for the problem becomes quite hairy:
% 
% $$2u_{i+1}-2u_{i}+u_{i}+u_{i-1}+hu_iu_{i+1}-hu_iu_{i-1}+2h^2u_i=0$$
%
% This problem is no longer linear, and cannot be solved with a sparse
% matrix solution.
%% Appendix
function f = newton2DR(fvec,x0,error,steps)
    syms x y
    v = [x,y];
    J = @(x) jacobian(fvec(x));
    Jeval = @(b) eval(subs(subs(J(v),b(1)),b(2)));
    
    if steps > 20
        f = 'error: poor inital guess or the solution does not converge.';
    elseif norm(fvec(x0))^2 <= error %stopping condition
        f = x0;
    else
        x0 = x0 - Jeval(x0)\fvec(x0);
        f = newton2DR(fvec,x0,error,steps+1);
    end
end
function f = newton2DRW(fvec,x0,error) % wrapper for recursion
    f = newton2DR(fvec,x0,error,0);
end

function f = newton2DI(fvec,x0,error) %iterave version of newton's method
    syms x y
    v = [x,y];
    J = @(x) jacobian(fvec(x));
    Jeval = @(b) eval(subs(subs(J(v),b(1)),b(2)));
    
    i=1;
    
while norm(fvec(x0))^2 > error
    x0 = x0 - Jeval(x0)\fvec(x0);
    i = i+1;
    if i > 20
        break;
    end
end
if i > 20
    f = 'error: poor guess, or solution does not converge';
else
    f = x0;
end
end

function f = newton2DMaxMinR(fvec,x0,error,steps)
    syms x y
    v = [x,y];
    
    gradF = @(x) gradient(fvec(x));
    gradFeval = @(b) eval(subs(subs(gradF(v),b(1)),b(2))); %derivatives need substitutions
    J = @(x) jacobian(gradF(x));
    Jeval = @(b) eval(subs(subs(J(v),b(1)),b(2)));
    
    if steps > 20
        f = 'error: poor inital guess or the solution does not converge.'
    elseif norm(gradFeval(x0))^2 <= error %stopping condition
        f = x0;
    else
        x0 = x0 - Jeval(x0)\gradFeval(x0);
        f = newton2DMaxMinR(fvec,x0,error,steps+1);
    end
end
function f = newton2DMaxMin(fvec,x0,error) % wrapper for recursion
    f = newton2DMaxMinR(fvec,x0,error,0);
end

function r = rungeKutta4(fun,y0,h)
    t = 0:h:2;
    y = zeros(size(t));
    y(1) = y0;
    [small big] = size(t);
    for i = 1:(big - 1)
        K1 = h * fun(t(i),y(i));
        K2 = h * fun(t(i)+h/2,y(i)+K1/2);
        K3 = h * fun(t(i)+h/2,y(i)+K2/2);
        K4 = h * fun(t(i)+h,y(i)+K3);
        y(i+1) = y(i) + (K1 + 2*K2 + 2*K3 + K4)/6;
    end
    r = y;
end