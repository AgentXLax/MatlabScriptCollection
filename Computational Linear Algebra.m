%% Assignment 4
% Nathan Wisla
%
% Numerical Methods (AUMAT 340)
% 
% March 14, 2016

%% 1) Continuous Least Squares (Midterm Re-Evaluation)
% Use the least-square method to fit a quadratic function to 
% fit $f(x) = sin{(\pi/2)}$ on $[0,2]$
% 
% Using the matrix of integrals from the lecture notes evaluated at [0,2],
% we obtain the following matrices to solve with in our linear system
% $Ax=b$.
%%
a = 0;
b = 2;
syms x;

g =@(x) sin(pi*x/2);

a11 = b - a;
a12 = int(x,[a b]);
a13 = int(x^2,[a b]);
a23 = int(x^3,[a b]);
a33 = int(x^4,[a b]);

A = [a11 a12 a13;a12 a13 a23;a13 a23 a33]

b1 = int(g(x), [a b]);
b2 = int(x*g(x), [a b]);
b3 = int(x^2*g(x), [a b]);
b = double([b1; b2; b3])

%%
% Our best-fitting coefficients and the resulting function give us values
% when the function is evaluated at the same points as it was in the
% midterm:

x = double(A\b); 
alpha = x(1)
beta = x(2)
gamma = x(3)
 
f =@(x) alpha + beta*x + gamma*x^2;
fAt0 = f(0)
fAt1 = f(1)
fAt2 =f(2)
%%
% we do not get the same result because the midterm problem was discrete
% and needed to fit exactly 3 points of $\sin{(\pi x/2)}$, where this problem is 
% continuous and needs to "best fit" *all* points of $\sin{(\pi x/2)}$.

%% 2) Non-Linear Least-Squares
% Use the least-square method to determine the constant of the function
% 
% $$g(x)=e^{bx}$$
% 
% so that it “best matches” the function 
% 
% $$f(x)=x^2 + 1$$
% 
% on $[0,1]$.
% 
% Using the least-squares method, we will arrive at the following equation:
% 
% $$0 = \int_{0}^{1}{(e^{bx} - (x^2 + 1))dx}$$
%%
a = 0
b = 1

syms c x;

f = @(x) x^2 + 1;
g = @(x) exp(c*x);

F = @(x) (g(x) - f(x))^2;
dFdc = @(x) diff(F(x),c);
%%
% Since we cannot analytically get our best value for $b$, we will put it
% in a graph and estimate the root.
%%
intEval =@(c) int(dFdc(x),[a b]);

fplot(intEval(x))
axis([-2 2 -1 5])
title('Values of b');
xlabel('b');
ylabel('y');

%%
% Fortunately for us, Matlab can already estimate the root of the function.
% We will plot it 'zoomed in' regardless to give a geometric interpretation
% of what's happening. Zooming in on the plot gives the following:
fplot(intEval(x))
axis([.604 .605 0 .001])
title('values of b (zoomed in)');
xlabel('b');
ylabel('y');
c = solve(intEval(x))

g =@(x) paramFunction(x,c);
%%
% We are left with the following plot, comparing the quadratic to our
% best-fit exponential:
%%

fplot(g(x),[0 1]);
hold on;
fplot(f(x),[0 1]);

title('x^2 + 1 vs. e^{bx}');
xlabel('x');
ylabel('y');

legend('x^2 + 1','e^{bx}');

hold off;

%% 3) L-U Decomposition
% Find the determinant of $A$, its inverse, and solve $Ax=b$ using L-U
% decomposition.
% 
% Using Nathan's algorithm (because he didn't look far enough into the
% lecture notes to find Doolittle's algorithm) We can obtain our lower and
% upper triangular matrices:
%%
A = [1 2 4;3 8 14;2 6 13]

L = lowerUpper(A,'lower') 
U = lowerUpper(A,'upper')

%%
% To confirm $LU = A$:
L*U

%%
% Using the fact that the determinant is a homomorphism under multiplication, 
% finding the determinant of an L-U decomposition only involves multiplying 
% the diagonals of the upper-matrix because $det(L)=1$:
%%
detLU = 1;
for i=1:3
    detLU = detLU*L(i,i)*U(i,i);
end    
detLU = detLU
%%
% To find the inverse of $A$ non trivially, I used the adjoint method of
% inverting both $L$ and $U$ by hand. Symbollically, it looks like the
% following:
%%
syms u11 u12 u13 u22 u23 u33 l21 l31 l32 detUsym

symInvL = [1 0 0;-l21 1 0;l21*l32-l31 -l32 1]
symInvU = 1/detUsym*[u22*u33 -u12*u33 u12*u23-u13*u22;0 u11*u33 -u11*u23;0 0 u11*u22]

%%
% I did this by hand, but I need to output my actual solutions in Matlab,
% so I'll do the lazy way now. The inverse of $A$ will be 
% $(LU)^{-1}=U^{-1}L^{-1}$
%%
InvL = inv(L)
InvU = inv(U)

inv(L*U)

%%
% Now to determine the system $Ly = b$ where $y = Ux$.
b = [1;2;3];

soln = U\(L\b)

x = rats(soln(1)) % rats() returns the answer as a fraction
y = rats(soln(2))
z = rats(soln(3))

%% 4) The Jacobi Method and the Gauss-Seidel Method
% Apply the Jacobi Method (and Gauss-Seidel Method for bonus marks) to
% solve
%
% $$5x-2y+3z=-1$$
% 
% $$-3x+9y+z=-1$$
% 
% $$2x-y-7z=-1$$
%
% Using the initial guess $x_0 = \langle0,0,0\rangle$
% and compare the methods with the exact solutions
%%
A = [5 -2 3;-3 9 1;2 -1 -7];
x0 = [0;0;0];
b = [-1;2;3];

%%
% these methods are in the function appendix section and solve for $x_n$
% recursively. An important thing to note is that I solved these values by
% using matrices instead of indices.
% 
% Jacobi Method: $\vec{x}_{n+1} = \hat{D}^{-1}(\vec{b}-(\hat{L}-\hat{U})\vec{x}_{n})$
% 
% Gauss-Seidel Method: $\vec{x}_{n+1} = \hat{L_{*}}^{-1}(\vec{b}-\hat{U})\vec{x}_{n})$
% 
% where $\hat{L_{*}}=\hat{L} + \hat{D}$

exactSolution = A\b
jacobiMethod = jacobi(A,b,x0,3) 
gaussSeidelMethod = gaussSeidel(A,b,x0,3)
%%
% We also obtain relative errors of:
%%
relJacobiError = abs(exactSolution - jacobiMethod)./exactSolution

relGaussError = abs(exactSolution - gaussSeidelMethod)./exactSolution

%% Function Appendix

function y = paramFunction(c,x)
    y = exp(c*x);
end % Q2 plugs in constant c to function handle

function LU = lowerUpper(aMatrix, lowerOrUpper);
    % i == j for diagonal
    % i > j for lower tri
    % i < j for upper tri
    U(1,1) = aMatrix(1,1);
    for j = 1:3;
        for i = 1:3;

            if i > j % below the diagonal
                U(i,j) = 0;
                sum = aMatrix(i,j);

                if j > 1 % not in the first column DO NOT TOUCH THE FIRST COLUMN
                    for k = 1:i-1;
                        sum = sum - L(i,k)*U(k,j);
                    end
                end
                L(i,j) = sum / U(j,j);

            elseif i <= j; % above and including the diagonal
                L(i,j) = 0;
                L(i,i) = 1;
                sum = aMatrix(i,j);

                for k = 1:i-1;
                    sum = sum - L(i,k)*U(k,j);
                end
                U(i,j) = sum; 

            end %if block
        end %row for
    end %column for
    if lowerOrUpper == 'lower';
        LU = L;
    elseif lowerOrUpper == 'upper';
        LU = U;
    end
    
end

function f = jacobi(A,b,x0,steps) % recursive function for jacobi's method
    for i = 1:3 % construct the diagonal inverse, lower and upper matrices
       for j = 1:3
           if i == j
              D(i,j) = 1/A(i,j); % inverse of a diagonal is 1/diagonal entries
              L(i,j) = 0;
              U(i,j) = 0;
           elseif i > j
               L(i,j) = A(i,j);
               D(i,j) = 0;
               U(i,j) = 0;
           else
               U(i,j) = A(i,j);
               L(i,j) = 0;
               D(i,j) = 0;
           end % conditionals
       end % column for
    end % row for

    % the meat of Jacobi's method done recursively
    if steps < 0
        f = 'error in number of steps (must be greater than or equal to 0)';
    elseif steps == 0
        f = x0;
    elseif steps == 1
        f = D*(b - (L + U)*x0);
    else 
        f = D*(b - (L + U)*jacobi(A,b,x0,steps - 1));
    end
end

function g = gaussSeidel(A,b,x0,steps)
    for i = 1:3
       for j = 1:3
           if i >= j
               L = A(i,j);
           else
               U = A(i,j);
           end
       end
    end

    if steps < 0
        g = 'step must be greater than zero'
    elseif steps == 0
        g = x;
    elseif steps == 1
        g = L\(b - U*x0); 
    else
        g = L\(b - U*gaussSeidel(A,b,x0,steps-1));
    end
end