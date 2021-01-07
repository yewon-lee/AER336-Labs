%% AER336: Assignment 6

% Yewon Lee
% 1004407260

%% Problem 2(c)

close all; clear; clc;

% define parameters
T = 50;
dt = 0.05;
u0 = [4;0];

% solve using RK4
odefun = @vdp;
[tt,uu] = rk4(odefun,T,dt,u0);

% Plot solution
figure(1)
plot(tt,uu(1,:));
xlabel('time');
ylabel('\phi(t)');
title('Q2(c): van der pol solution using RK4')
hold off

%% Problem 2(d)

close; clear all; clc;

% define parameters
T = 50;
dt = 0.05;
u0 = [4;0];

% solve using RK4
odefun = @vdp;
[tt,uu] = rk4(odefun,T,dt,u0);

% compute eigenvalues
for j = 1:length(tt)
    [~,dfdu] = vdp(tt(j),uu(:,j));
    lambda(2*j-1:2*j) = eig(dfdu);
end

% plot eigenvalues
figure(2)
plot(real(lambda),imag(lambda),'x')
xlabel('Real')
ylabel('Imaginary')
title('Q2(d): eigenvalues of dfdu')
hold off

%% Problem 2(e)

close all; clear; clc;

% parameters
T = 50;
u0 = [4;0];
dt1 = 0.1;
dt2 = 0.05;
dt3 = 0.025;
dt4 = 0.0125;

% reference solution (Dormand Prince)
odefun = @vdp;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[tref,uref] = ode45(odefun,[0 T],u0,opts);
tref = tref';
uref = uref';

% approximate solutions (RK4)
[tt1,uu1] = rk4(odefun,T,dt1,u0);
[tt2,uu2] = rk4(odefun,T,dt2,u0);
[tt3,uu3] = rk4(odefun,T,dt3,u0);
[tt4,uu4] = rk4(odefun,T,dt4,u0);

% compute error
e1 = abs(uu1(1,end) - uref(1,end))
e2 = abs(uu2(1,end) - uref(1,end))
e3 = abs(uu3(1,end) - uref(1,end))
e4 = abs(uu4(1,end) - uref(1,end))

% plot error and calculate convergence rate
figure(3)
e = [e1 e2 e3 e4];
dts = [dt1 dt2 dt3 dt4];
loglog(dts,e,'-o');
title('Q2(e): loglog plot of error versus time step');
p = polyfit(log(dts(2:end)), log(e(2:end)), 1);
conv_rate_rk4 = p(1)

% COMMENTS: 
% The observed convergence rate is 3.9494 which is ~4. This is
% consistent with our expectations since RK4 is theoretically fourth order 
% convergent.
% The behavior we observe with dt = 0.1 is instability (our numerical solution
% "blows" up over time). As we can see in the real-imaginary plot of
% eigenvalues of df/du at different points in time, the eigenvalues can
% grow very large for certain t (for example, one of the eigenvalues is at
% about -45). For such eigenvalues, lambda*dt end up outside of the RK4
% stability region (for e.g. -45*0.1 = -4.5, which is indeed outside of the RK4
% stability region). This will cause the scheme to become unstable.

%% Problem 3(b)

clear; close all; clc;

% Compute matrix B
n1 = 16;
n = n1^2;
A = laplacian2d(n1);
B = [sparse(n,n), speye(n); A, sparse(n,n)];

% Compute eigenvalues of B
lambda = eig(full(B));

% Plot eigenvalues on complex plane
figure(4)
plot(real(lambda),imag(lambda),'o')
xlabel('Real')
ylabel('Imaginary')
title('Q3(b): eigenvalues of B')
daspect([1 1 1])
pbaspect([1 1 1])
hold off

% COMMENTS:The eigenvalues of B all lie on the imaginary axis (i.e. the
% eigenvalues are all imaginary). This is expected; the wave equation
% admits traveling waves OR standing waves (in our case, standing waves 
% because of the homogeneous Dirichlet BC) as a solution. These waves do
% not disspiate energy or decay/grow. Therefore, we do not expect a positive
% real component of lambda (implying growing amplitude) or negative real
% component of lambda (decaying amplitude).

%% Problem 3(c)

clear; close all; clc;

% Compute matrix B
n1=128; 
n = 16384;
A = laplacian2d(n1);
B = [sparse(n,n), speye(n); A, sparse(n,n)];

% Compute largest eigenvalue
big_eigs = eigs(B,6,'largestimag','Tolerance',1e-4);
max_eig = max(big_eigs)
magnitude_max_eig = abs(max_eig)

% COMMENTS: the magnitude of the maximum eigenvalue is 364.7045

%% Problem 3(d)

clear; close all; clc;

% Parameters
n1 = 128;
n = n1^2;
odefun = @wave;
dt = 0.005; 
T = 0.3;
load('prob3ics.mat')
u0 = [x0;v0];

% Solve wave equation using RK4
[tt,uu] = rk4_mod(odefun,T,dt,u0);

% Plot solution at different times
figure(5)
plotfield2d(uu(1:n,1));
title('Q3(d): solution at t = 0')
xlabel('x');
ylabel('y')
hold off

figure(6)
plotfield2d(uu(1:n,0.1/dt+1));
%tt(0.1/dt+1)
title('Q3(d): solution at t = 0.1')
xlabel('x');
ylabel('y')
hold off

figure(7)
plotfield2d(uu(1:n,0.2/dt+1));
%tt(0.2/dt+1)
title('Q3(d): solution at t = 0.2')
xlabel('x');
ylabel('y')
hold off

figure(8)
plotfield2d(uu(1:n,0.3/dt+1));
%tt(0.3/dt+1)
title('Q3(d): solution at t = 0.3')
xlabel('x');
ylabel('y')
hold off

% JUSTIFICATION (choice of dt):
% From the stability region of RK4, we know that for eigenvalues on the 
% imaginary axis, abs(lambda*dt) < 3 (approximately) to ensure that we are
% on/within the staiblity boundary. Because the actual limit is less than
% 3, I decided to use abs(lambda*dt) < 2.5 to choose dt, were lambda is 
% the maximum magnitude of all eigenvalues of B (364.7045 in our case). From
% abs(lambda)*dt=2.5, we get dt = 0.0069. I chose dt = 0.005 so that t =
% 0.1, 0.2, 0.3 could be obtained through time marching.

%% Problem 3(e)

clear; close all; clc;

% Parameters
n1 = 128;
n = n1^2;
odefun = @wave;
dt = 0.005; 
T = 0.3;
load('prob3ics.mat')
u0 = [x0;v0];

% Solve wave equation using RK4
[tt,uu] = rk4_mod(odefun,T,dt,u0);

% Plot at t=0.3
figure(9)
plotfield2d(uu(1:n,0.3/dt+1));
title('Q3(d): solution at t = 0.3')
xlabel('x');
ylabel('y')
view(0,90)
hold off

% COMMENT: the six letter message is 'AER336' :)

%% Problem 2(a)

function [f,dfdu] = vdp(t,u)
    mu = 3;
    u1 = u(1);
    u2 = u(2);
    f = [u2 ; mu*(1-u1^2)*u2 - u1 + sin(t)];
    dfdu = [0 , 1; -2*mu*u1*u2 - 1, mu*(1-u1^2)];
end

%% Problem 2(b)

function [tt,uu] = rk4(odefun,T,dt,u0)
    J = T/dt + 1; % number of time steps plus one for IC
    uu = NaN(2,J);
    tt = (0:dt:T)';
    uu(:,1) = u0;
    for j = 2:J
        v1 = uu(:,j-1);
        [F1,~] = odefun(tt(j-1), v1);
        v2 = uu(:,j-1) + 0.5*dt*F1;
        [F2,~] = odefun(tt(j-1) + 0.5*dt, v2);
        v3 = uu(:,j-1) + 0.5*dt*F2;
        [F3,~] = odefun(tt(j-1) + 0.5*dt, v3);
        v4 = uu(:,j-1) + dt*F3;
        [F4,~] = odefun(tt(j-1) + dt, v4);
        uu(:,j) = uu(:,j-1) + dt*(F1/6 + F2/3 + F3/3 + F4/6);
    end
end


%% Other function definitions

% Q3: 2d Laplacian centered difference matrix A
function A = laplacian2d(n1) 
    h = 1/(n1+1); 
    A1 = 1/h^2*spdiags([ones(n1,1),-2*ones(n1,1),ones(n1,1)],[-1,0,1],n1,n1); 
    E1 = speye(n1); 
    A = kron(A1,E1) + kron(E1,A1); 
end

% Q3(d): odefunction for wave equation
function [f,dfdu] = wave(t,u)
    % Parameters
    n1 = 128;
    n = n1^2;
    A = laplacian2d(n1);
    B = [sparse(n,n), speye(n); A, sparse(n,n)];
    
    % compute f and dfdu
    f = B*u;
    dfdu = 0; % we won't need this anyway
end

% Q3(d): modified RK4
function [tt,uu] = rk4_mod(odefun,T,dt,u0) 
    J = T/dt + 1; % number of time steps plus one for IC
    n = length(u0);
    uu = NaN(n,J);
    tt = (0:dt:T)';
    uu(:,1) = u0;
    for j = 2:J
        v1 = uu(:,j-1);
        [F1,~] = odefun(tt(j-1), v1);
        v2 = uu(:,j-1) + 0.5*dt*F1;
        [F2,~] = odefun(tt(j-1) + 0.5*dt, v2);
        v3 = uu(:,j-1) + 0.5*dt*F2;
        [F3,~] = odefun(tt(j-1) + 0.5*dt, v3);
        v4 = uu(:,j-1) + dt*F3;
        [F4,~] = odefun(tt(j-1) + dt, v4);
        uu(:,j) = uu(:,j-1) + dt*(F1/6 + F2/3 + F3/3 + F4/6);
    end
end



