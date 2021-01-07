%% Assignment 4

% Yewon Lee
% 1004407260

%% Problem 2 (b)

close; clear all; clc;

y = @fun;

% test the bisection function
a0 = 2; b0 = 0;
tol = 10^-4;
[a,b] = my_bisection(y,a0,b0,tol)

% Comments: We get a = 1.0001 and b = 1, which is the correct root of our
% funciton to within the tolerance of 10^-4.


%% Problem 2 (c) Part 2: Confirm function defined in Part 1 works

close; clear all; clc;

% test the bisection function
x0 = 5;
tol = 10^-4;
y = @fun;
x = my_newton(y,x0,tol) % x value

% Comment: Testing the function, we get x = 1.0000 to machine precision
% which is correct.

%% Problem 2 (d)

g = @(x) nonlinear_bench(x);

% Bisection method
a0 = 0; b0 = 5; tol = 10^-4;
[a_bench,b_bench] = my_bisection(g,a0,b0,tol)

% Newton's method
x0 = 0;
x_bench = my_newton(g,x0,tol)

%% Problem 2 (e)

g = @(x) nonlinear_bench(x);

% Newton's method
x0 = i; tol = 10^-4;
x_bench = my_newton(g,x0,tol)

%% Problem 3 (c)

clear; close all; clc;

x0 = [0;0];
tol = 10^-4;
y = @ my_fun;
[x,fmin] = my_newton_opt(y,x0,tol)

% Comment: The minimum occurs at (1,1) and the value of the function at
% that point is 0.


%% Problem 3 (d) Part 2

clear; close all; clc;

x0 = [0.2;0.2];
tol = 10^-4;
y = @ my_fun;
[x,xi,fmin] = my_newton_opt_mod(y,x0,tol)
err = NaN(length(xi(1,:)),1);
for i = 1:length(err)
    err(i) = norm(x-xi(:,i));
    %err(i) = (x-xi(:,i))' * (x-xi(:,i));
end

err % each ith index corresponds to time i-1, with i=1,2,3,...

% Comments: 
% Except for the first ~3 entries (note that the first entry corresponds to
% the initial guess so can be ignored), the error converges pretty much
% quadratically. For example, for err(i=5) = 0.7081, and 0.7081^2 = 0.501
% which is close to err(i=6) = 0.3822. Similarly, err(i=7) is close to 
% 0.3822^2 = 0.1461.

%% Problem 3 (e)

clear; close all; clc;

% Newton's method
x0 = [0;1.0];
tol = 10^-4;
y = @ my_fun;
[x,xi,fmin] = my_newton_opt_mod(y,x0,tol)

% Contour plot
x = linspace(-1.25,1.25);
y = linspace(-0.75,1.75);
[X,Y] = meshgrid(x,y);
Z = (1-X).^2 + 10*(Y-X.^2).^2;
contourf(X,Y,Z,50)

hold on
plot(xi(1,:),xi(2,:),'o-r')

%% Problem 2 (a)

function [a,b] = my_bisection(fun,a0,b0,tol)
    a = a0; b = b0;
    err = abs(b-a);
    while err > tol
        c = (a + b)/2;
        [fa, dfa] = fun(a);
        [fb, dfb] = fun(b);
        [fc, dfc] = fun(c);
        if sign(fc) == sign(fa)
            a = c;
        else
            b = c;
        end
        err = abs(b-a);
    end
end

%% Problem 2 (c) Part 1

function x = my_newton(fun,x0,tol)
    [f, df] = fun(x0);
    err = abs(f);
    x = x0;
    [f, df] = fun(x);
    while err > tol
        x = x - f/df;
        [f, df] = fun(x);
        err = abs(f);
    end
end

%% test function for Problem 2

function [f,df] = fun(x)
    y = @(x) (x-1)*(x^2 + 1); % the real root is x=1
    y_prime = @(x) 3*x^2 - 2*x + 1;
    f = y(x);
    df = y_prime(x);
end

%% Problem 3 (a)

function [f,df,ddf] = my_fun(z)
    x = z(1); y = z(2);
    f = (1-x)^2 + 10*(y-x^2)^2;
    
    df = [2*(x-1) - 40*x*(y-x^2); 20*(y-x^2)];
      
    ddf = [2-40*y+120*x^2, -40*x; -40*x, 20];
end

%% Problem 3 (b)

function [x,fmin] = my_newton_opt(my_fun,x0,tol)
    x = x0;
    [f,df,ddf] = my_fun(x);
    err = tol + 1; % Just set initial error to be bigger than tol so that we always enter the while loop
    while err > tol
        dx = -ddf\df;
        x = x + dx;
        err = dx' * dx;
        [f,df,ddf] = my_fun(x);
    end
    fmin = f;
end

%% Problem 3 (d) Part 1

% Modification of my_newton_opt

function [x,xi,fmin] = my_newton_opt_mod(my_fun,x0,tol)
    x = x0;
    xi = [x];
    [f,df,ddf] = my_fun(x);
    err = tol + 1; % Just set initial error to be bigger than tol so that we always enter the while loop
    while err > tol
        dx = -ddf\df;
        x = x + dx;
        xi = [xi,x];
        err = dx' * dx;
        [f,df,ddf] = my_fun(x);
    end
    fmin = f;
end

