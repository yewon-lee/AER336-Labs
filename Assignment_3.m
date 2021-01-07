%% AER 336 Assignment 3

% Yewon Lee
% 1004407260

%% Problem 2 (b)

close all; clear; clc;

% test function myforwardsub

% L matrix 
L = [3 0 0;
    1 1 0;
    2 1 1];

% b vector
b = [3 -1 1]';

% solve for x
x = myforwardsub(L,b)

%% Problem 2 (e)

close all; clear; clc;


% A matrix
A = [2 1 0;
    -1 2 -1;
    0 -1 3];

% perform LU factorization
[L,U] = mylu(A)

%% Problem 3 (a)

close all; clear; clc;


% Parameters
m3 = 3; m5 = 5; m9 = 9;

% Generate matrices
[A3, b3] = poisson.getmatvec(m3,1);
[A5, b5] = poisson.getmatvec(m5,1);
[A9, b9] = poisson.getmatvec(m9,1);

format long

% check the eigenvalues of the A matrices
d3 = eig(A3)
d5 = eig(A5)
d9 = eig(A9)

% check the eigenvalues of the A matrices with rounding
d3_rounded = round(eig(A3))
d5_rounded = round(eig(A5))
d9_rounded = round(eig(A9))

% COMMENT ON SPD: based on the UNROUNDED eigenvalues of A3, A5, and A9, A3 and A5 are
% SPD (which correspond to when m=3 and m=5). A9 (which corresponds to m=9)
% is not SPD, because it has complex eigenvalues when UNROUNDED.
% However, with rounding, all three matrices are SPD, as they have positive
% eigenvalues. It appears that the imaginary parts of A9's eigenvalues are
% relatively small (~machine precision) compared to the real part. Thus,
% all three matrices are pretty much SPD.

%% Problem 3 (b)

close all; clear; clc;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  PART 1: m=3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% case m=3, flag=0
m = 3; flag = 0;
[A0,b0] = poisson.getmatvec(m,flag);
x0 = A0\b0;
figure(1)
poisson.vizsoln(x0,flag)
title('m=3, flag=0')
hold off

% case m=3, flag=1
flag = 1;
[A1,b1] = poisson.getmatvec(m,flag);
x1 = A1\b1;
figure(2)
poisson.vizsoln(x1,flag)
title('m=3, flag=1')
hold off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  PART 1: m=9  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% case flag=0
m = 9; flag = 0;
[A0,b0] = poisson.getmatvec(m,flag);
x0 = A0\b0;
figure(3)
poisson.vizsoln(x0,flag)
title('m=9, flag=0')
hold off

% case flag=1
flag = 1;
[A1,b1] = poisson.getmatvec(m,flag);
x1 = A1\b1;
figure(4)
poisson.vizsoln(x1,flag)
title('m=9, flag=1')
hold off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  PART 1: m=18  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% case flag=0
m = 18; flag = 0;
[A0,b0] = poisson.getmatvec(m,flag);
%x0 = A0\b0;
x0 = mylinearsolver(A0,b0);
figure(5)
poisson.vizsoln(x0,flag)
title('m=18, flag=0')
hold off

% case flag=1
flag = 1;
[A1,b1] = poisson.getmatvec(m,flag);
x1 = A1\b1;
%x1 = mylinearsolver(A1,b1);
figure(6)
poisson.vizsoln(x1,flag)
title('m=18, flag=1')
hold off


% COMMENTS
% Comparing Figures 1 and 2, we see that the solution to m=3 and flag=0/1
% are indeed identical. The same goes for the solution to m=9 and flag=0/1
% shown in Figures 3 and 4.
% However, the solution to m=18 and flag=0/1 are visibly different. In
% fact,for flag=0, we incur an error message saying that the matrix is close to
% singular with RCOND = 2.474792e-20 (i.e. the matrix is ill-conditioned).
% The ill conditioning of this matrix is probably the reason for this
% discrepancy. 

%% Problem 3 (c)

close all; clear; clc;

% assess the computational cost of mylu algorithm

flag = 1; % stable flag
m_sq = (10:2:30)'; 
m_sq = m_sq.^2;
t = NaN(length(m_sq),1);

for m = 10:2:30
    [A,b] = poisson.getmatvec(m,flag);
    tic;
    [L U] = mylu(A);
    t(int8(1+(m-10)/2)) = toc;
end

figure(7)
loglog(m_sq,t);
xlabel('m^2'); ylabel('t');
title('Computational Cost vs. m^2 of LU Factorization');
p = polyfit(log(m_sq(2:end)), log(t(2:end)), 1);
slope_mylu = p(1)
hold off

% COMMENTS: 
% We know from our analysis that the cost of LU factorization scales
% approximately by (2/3)n^3 for large n, where n is m^2 in our case. 
% Upon plotting computation time as a function of matrix size (m^2), I got 
% the slope of the log-log plot to be 2.8428, which is about 3. 3 is the 
% power on n that we expect, so our observations are somewhat in line with
% theoretical predictions. 
% Note that the slope changes every time I run the code, so the printed
% slope might differ from what I report above.

%% Problem 3 (d)

close all; clear; clc;

% assess computational cost of myforwardsub algorithm

flag = 1; % stable flag
m = (10:2:30)'; 
m_sq = m.^2;
t = NaN(length(m_sq),1);

for i = 1:length(m)
    [A,b] = poisson.getmatvec(m(i),flag);
    [L,U] = mylu(A);
    tic;
    x = myforwardsub(L,b);
    t(i) = toc;
end

figure(8)
loglog(m_sq(2:end),t(2:end),'-o');
xlabel('m^2'); ylabel('t');
title('Computational Cost vs. m^2 of Forward Substitution');
p = polyfit(log(m_sq(3:end)), log(t(3:end)), 1);
slope_myforwardsub = p(1)
hold off

% COMMENTS: 
% We know that the relationship between cost and matrix size n follows approximately
% n^2 for forward substitution. In our case, n is equal to m^2. Plotting
% the logarithmic size of the matrix (m^2) against the logarithmic cost, I
% got the slope to be 1.46055. The expected value of the slope is ~2, so
% this is approximately in line with theoretical expectations.


%% Problem 4 (a),(b),(c),(d)

close all; clear; clc;

% Parameters
dt = 0.05;
T = 10;
w = 2;

% assemble matrices and vectors
Y = signal_sampler(dt,T);
t = (0:dt:T)';

m = length(Y); n = 3;
h0 = @(t) cos(w*t);
h1 = @(t) sin(w*t);
X = ones(m,n);
for i = 1:m
    X(i,1) = h0(t(i));
    X(i,2) = h1(t(i));  
end

% Compute max likelihood estimates
beta = (X' * X)\(X' * Y);

% Find a,b,c
a = beta(1)
b = beta(2)
c = beta(3)

% COMMENTS:
% The assumptions I made on the noise is that it is (i) normal with zero
% mean, (ii) homoscedastic or that the distribution is not dependent on x,
% and (ii) the noise is independent. I also assumed that the model is
% unbiased.

% Plot regression
Ymle = @(t) a*cos(w*t) + b*sin(w*t) + c;
figure(9)
ezplot(Ymle,[0,10]);
hold on
plot(t,Y,'o');
xlabel('t');
ylabel('Y');
legend('Y_{model}','Y')
hold off

% Estimate std dev of measurement noise
std_dev = (X*beta - Y)' * (X*beta - Y) / (m-n);
std_dev = sqrt(std_dev)

% 95% confidence interval for a,b,c
sigma = (std_dev^2) * inv(X' * X);
t95 = tinv(0.5 + 0.95/2,m-n);
Ia = [beta(1) - t95*sqrt(sigma(1,1)), beta(1) + t95*sqrt(sigma(1,1))]
Ib = [beta(2) - t95*sqrt(sigma(2,2)), beta(2) + t95*sqrt(sigma(2,2))]
Ic = [beta(3) - t95*sqrt(sigma(3,3)), beta(3) + t95*sqrt(sigma(3,3))]

% half width of each confidence interval
hw_a = t95*sqrt(sigma(1,1))
hw_b = t95*sqrt(sigma(2,2))
hw_c = t95*sqrt(sigma(3,3))

%% Problem 4 (e)

close all; clear; clc;

% COMMENTS:
% The confidence interval for each parameter scales roughly by a factor of
% 1/sqrt(m-n). Since n=3 in our case and our m is much larger at 201, this 
% scale can again be approximated by ~1/sqrt(m). Hence, to decrease the
% confidence interval by a factor of 10, we need to increase m by a factor
% of 100 (that is, use approx 20100 data points). With dt fixed, this means
% that our T needs to be 100 times longer than before.

% Parameters
dt = 0.05;
T = 1000;
w = 2;

% assemble matrices and vectors
Y = signal_sampler(dt,T);
t = (0:dt:T)';

m = length(Y); n = 3;
h0 = @(t) cos(w*t);
h1 = @(t) sin(w*t);
X = ones(m,n);
for i = 1:m
    X(i,1) = h0(t(i));
    X(i,2) = h1(t(i));  
end

% Compute max likelihood estimates
beta = (X' * X)\(X' * Y);

% Find a,b,c
a = beta(1)
b = beta(2)
c = beta(3)

% Estimate std dev of measurement noise
std_dev = (X*beta - Y)' * (X*beta - Y) / (m-n);
std_dev = sqrt(std_dev);

% Recompute the 95% confidence interval for a,b,c
sigma = (std_dev^2) * inv(X' * X);
t95 = tinv(0.5 + 0.95/2,m-n);
Ia = [beta(1) - t95*sqrt(sigma(1,1)), beta(1) + t95*sqrt(sigma(1,1))]
Ib = [beta(2) - t95*sqrt(sigma(2,2)), beta(2) + t95*sqrt(sigma(2,2))]
Ic = [beta(3) - t95*sqrt(sigma(3,3)), beta(3) + t95*sqrt(sigma(3,3))]

% half width of each confidence interval
hw_a = t95*sqrt(sigma(1,1))
hw_b = t95*sqrt(sigma(2,2))
hw_c = t95*sqrt(sigma(3,3))

% COMMENTS:
% In the previous section, our halfwidths were 
% hw_a = 0.0378
% hw_b = 0.0385
% hw_c = 0.0270
%
% With T = 1000, we get the halfwidths to be
% hw_a = 0.0039
% hw_b = 0.0039
% hw_c = 0.0028
%
% Indeed, the confidence intervals were reduced by approximately a factor
% of 10 by setting T=1000.



%% Problem 2 (a)

function x = myforwardsub(L,b)
    n = length(b);
    x = NaN(n,1);
    x(1) = b(1)/L(1,1);
    for i = 2:n
        x(i) = b(i);
        for j = 1:i-1
            x(i) = x(i) - L(i,j)*x(j);
        end
        x(i) = x(i) / L(i,i);
    end
end

%% Problem 2 (c)

function x = mybackwardsub(U,b)
    n = length(b);
    x = zeros(n,1);
    x(n) = b(n)/U(n,n);
    for i = n-1:-1:1
        x(i) = b(i);
        for j = i+1:n
            x(i) = x(i) - U(i,j)*x(j);
        end
        x(i) = x(i)/U(i,i);
    end
end


%% Problem 2 (d)

function [L,U] = mylu(A)
    U = A;
    n = length(A(:,1));
    L = eye(n);
    for i = 1:n-1
        for j = i+1:n
            L(j,i) = U(j,i)/U(i,i);
            for k = i:n
                U(j,k) = U(j,k) - L(j,i)*U(i,k);
            end
        end
    end
end
    
%% Problem 2 (f)

function x = mylinearsolver(A,b)
    [L,U] = mylu(A); 
    z = myforwardsub(L,b);
    x = mybackwardsub(U,z);
end

