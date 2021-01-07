%% AER336 Assignment 2

% Yewon Lee
% 1004407260

%% Problem 2 (b)

% Functions to integrate
f0 = @(x) 1;
f1 = @(x) x;
f2 = @(x) x^2;
f3 = @(x) x^3;
f4 = @(x) x^4;
f5 = @(x) x^5;

% Analytical solutions of integrals
I = [1 1/2 1/3 1/4 1/5 1/6]';

% Parameter definitions
a = 0; b = 1;

% Evaluate functions with composite_gauss2

% n = 1
n = 1;
q10 = composite_gauss2(f0,a,b,n);
q11 = composite_gauss2(f1,a,b,n);
q12 = composite_gauss2(f2,a,b,n);
q13 = composite_gauss2(f3,a,b,n);
q14 = composite_gauss2(f4,a,b,n);
q15 = composite_gauss2(f5,a,b,n);
q1 = abs([q10; q11; q12; q13; q14; q15] - I);

% n = 2
n = 2;
q20 = composite_gauss2(f0,a,b,n);
q21 = composite_gauss2(f1,a,b,n);
q22 = composite_gauss2(f2,a,b,n);
q23 = composite_gauss2(f3,a,b,n);
q24 = composite_gauss2(f4,a,b,n);
q25 = composite_gauss2(f5,a,b,n);
q2 = abs([q20; q21; q22; q23; q24; q25] - I);

% n = 4;
n = 4;
q40 = composite_gauss2(f0,a,b,n);
q41 = composite_gauss2(f1,a,b,n);
q42 = composite_gauss2(f2,a,b,n);
q43 = composite_gauss2(f3,a,b,n);
q44 = composite_gauss2(f4,a,b,n);
q45 = composite_gauss2(f5,a,b,n);
q4 = abs([q40; q41; q42; q43; q44; q45] - I);

% n = 8;
n = 8;
q80 = composite_gauss2(f0,a,b,n);
q81 = composite_gauss2(f1,a,b,n);
q82 = composite_gauss2(f2,a,b,n);
q83 = composite_gauss2(f3,a,b,n);
q84 = composite_gauss2(f4,a,b,n);
q85 = composite_gauss2(f5,a,b,n);
q8 = abs([q80; q81; q82; q83; q84; q85] - I);

% Construct a table of the error values
n_vec = [1 2 4 8]';
table_matrix = [q1 q2 q4 q8]';
table_matrix = [n_vec table_matrix];
Table = array2table(table_matrix, 'VariableNames', {'n', 'f0', 'f1', 'f2', 'f3', 'f4','f5'})

% COMMENT on observed error behavior:
%
% The Gauss quadrature rule integrates exactly polynomials of order up to
% 2m-1, where m is the number of points. In our case, since m=2, this would
% mean that we could integrate at most polynomials of degree 3 exactly.
% Observing the error results in the table below, we can see that this is
% indeed the case, since the errrors are zero (or machine zero) up until 
% f3 = x^3. 
%
% We observe that for polynomials of degree 4 and greater (i.e. those that 
% cannot be integrated exactly), the convergence rate is 4. That is,
% every 2-fold reduction in segment size results in a 16-fold reduction
% in the error.


%% Problem 3 (b)

% Define function to integrate and parameters
a = -1; b = 1;
f = @(x) 2*sqrt(1 - x^2);
I = pi; % actual value of the integral

% tol = 10^-2
tol2 = 10^-2;
[q2,e2,cnt2] = adaptive_gauss_kronrod(f,a,b,tol2);
e_actual2 = abs(I - q2);

% tol = 10^-5
tol5 = 10^-5;
[q5,e5,cnt5] = adaptive_gauss_kronrod(f,a,b,tol5);
e_actual5 = abs(I - q5);

% tol = 10^-8
tol8 = 10^-8;
[q8,e8,cnt8] = adaptive_gauss_kronrod(f,a,b,tol8);
e_actual8 = abs(I - q8);

% Construct a table of the error values
tols = [tol2 tol5 tol8]';
Q = [q2 q5 q8]';
error_actual = [e_actual2 e_actual5 e_actual8]';
error_estimate = [e2 e5 e8]';
cnt = [cnt2 cnt5 cnt8]';
table_matrix = [tols Q error_actual error_estimate cnt];
Table = array2table(table_matrix, 'VariableNames', {'tols', 'Q', 'error_actual', 'error_estimate', 'cnt'})

% COMMENT on the accuracy of the error estimate:
% As seen in the table, the error estimate is quite pessimistic; it is
% greater than the actual error by an entire order of magnitude, and even
% more, in the case of tols = 0.01. This makes sense, because the error
% estimate is an upper bound for the error.

%% Problem 3 (c)

% Compute the volume using adaptive_gauss_kronrod
format long

% NOTE: I will assume that the question asks for 4 sig fig accuracy WITH 
% ROUNDING (i.e. correct volume to four sigfigs after the computed integral
% is rounded)

tol = 10^-6; a = 0; b = 1; 
y = @ktgeom;
[q,e,cnt] = adaptive_gauss_kronrod(y,a,b,tol);
volume_AGK = 2*q % NOTE: AGK = adaptive gauss kronrod
error_in_integral_AGK = e % error estimate from using AGK
error_in_volume_AGK = 2*e % error estimate in the entire volume is twice the error in the integral approximation
func_queries = cnt % number of function queries using AGK

% Compute the volume using composite_gauss2
n = 1000;
q = composite_gauss2(y,a,b,n);
volume_composite_gauss2 = 2*q


%% Problem 2 (a)

% use the 2 point Gauss rule over each of the n segments

function q = composite_gauss2(f,a,b,n)

% 2 point Gauss quadrature weights and points on [-1,1]
x1 = -1/sqrt(3); x2 = 1/sqrt(3);
w1 = 1; w2 = 1;

% Initialize vector of x and w values over [a,b]
x = NaN(2*n,1); w = NaN(2*n,1);

% Map from [-1,1] to [a_new, b_new] 

for i = 1:n % loop over segments
    a_new = (b-a)/n * (i-1) + a;
    b_new = a_new + (b-a)/n;
    x(2*(i-1) + 1) = a_new + (b_new - a_new)/2 * (x1 + 1);
    x(2*i) = a_new + (b_new - a_new)/2 * (x2 + 1);
    w(2*(i-1) + 1) = (b_new - a_new)/2 * w1;
    w(2*i) = (b_new - a_new)/2 * w2;
end

F = NaN(2*n,1);
for i = 1:2*n
    F(i) = f(x(i));
end

q = F' * w;

end

%% Problem 3 (a)

% Adaptive quadrature

function [q,e,cnt] = adaptive_gauss_kronrod(f,a,b,tol)

% Load the quadrature points and weights
load('gauss_kronrod.dat');
x = gauss_kronrod(:,1); % quadrature points on [-1,1]
w_gauss = gauss_kronrod(:,2); % 7-point Gauss weights
w_kronrod = gauss_kronrod(:,3); % 15-point Kronrod weights

c = (a + b)/2;

% Compute Q0 (using 7-point Gauss quadrature)
Q0 = 0;
for i = 1:7 % loop over gauss quad points
    x_new = a + (b-a)/2 * (x(2*i) + 1);
    w_new = (b-a)/2 * w_gauss(2*i);
    Q0 = Q0 + f(x_new) * w_new;
end

% Compute Q1 (using 15-point Kronrod weights)
Q1 = 0;
for i = 1:15 % loop over Kronrod quad points
    x_new = a + (b-a)/2 * (x(i) + 1);
    w_new = (b-a)/2 * w_kronrod(i);
    Q1 = Q1 + f(x_new) * w_new;
end

% define cnt
cnt = 15;
    
% define the error, e
e = (200 * abs(Q1 - Q0))^(3/2);

% error tolerance condition
if e <= tol
    q = Q1;
else
    [qL,eL,cntL] = adaptive_gauss_kronrod(f,a,c,tol/2);
    [qR,eR,cntR] = adaptive_gauss_kronrod(f,c,b,tol/2);
    q = qL + qR;
    e = eL + eR;
    cnt = cntL + cntR - 1;
end
    
end



