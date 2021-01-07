%% AER336 Assignment 5

% Yewon Lee
% 1004407260

%% Problem 2(d)

close all; clear; clc;

% parameters
T = 10; 
dt1 = 1;
dt2 = 0.5;
dt3 = 0.25;

% solutions using forward Euler
func = @odefun;
[tt1,uu1] = forward_euler(func,T,dt1);
[tt2,uu2] = forward_euler(func,T,dt2);
[tt3,uu3] = forward_euler(func,T,dt3);

% plot solutions
figure(1)
plot(tt1,uu1);
hold on
plot(tt2,uu2);
plot(tt3,uu3);
syms t
ezplot(cos(t),[0,10]);
xlabel('time'); 
ylabel('u');
title('Q2(d): Forward Euler solutions');
legend('dt=1','dt=0.5','dt=0.25','analytical')
hold off

%% Problem 2(e)

close all; clc; clear;

% time steps and parameters
T = 10; 
dt1 = 1;
dt2 = 0.5;
dt3 = 0.25;
dt4 = 0.125;
dt5 = 0.0625;

% compute error at the final time
func = @odefun;
[tt1,uu1] = forward_euler(func,T,dt1);
[tt2,uu2] = forward_euler(func,T,dt2);
[tt3,uu3] = forward_euler(func,T,dt3);
[tt4,uu4] = forward_euler(func,T,dt4);
[tt5,uu5] = forward_euler(func,T,dt5);
err1 = abs(cos(10) - uu1(end));
err2 = abs(cos(10) - uu2(end));
err3 = abs(cos(10) - uu3(end));
err4 = abs(cos(10) - uu4(end));
err5 = abs(cos(10) - uu5(end));

% log-log plot of error vs. timestep
dts = [1, 0.5, 0.25, 0.125,0.0625]';
errs = [err1, err2, err3, err4, err5]';
figure(2)
loglog(dts,errs,'-o','LineWidth',1.5);
title('Q2(e): Convergence rate of Forward Euler');
xlabel('dt');
ylabel('error at T=10');
hold off

% convergence rate of FE
p = polyfit(log(dts),log(errs),1);
conv_rate_FE = p(1)

% COMMENTS: We get the the convergence rate of the forward euler method is
% approximately 1. This is consistent with expectations, since the FE
% method is first order accurate.

%% Problem 3(b)

close all; clc; clear;

% time steps and parameters
T = 1; 
dt = 0.5;

% compute error at the final time
func = @odefun_3b;
[tt,uu] = crank_nicolson(func,T,dt);
err_CN = abs(2 - uu(end)) % 2 is the analytical solution at time T=1

% COMMENTS: The error of the Crank-Nicolson method is bounded by a product
% consisting of the third derivative of u. In our case, the third
% derivative of u is 0, so the error should be zero. Experimenting with various
% time steps, we see that the error is indeed always "zero" (for any dt) 
% in the sense that the error is the same order of magnitude (or less) than
% the tolerance set by Newton's Method. In our case, the error is very close
% to machine precision. This is true regardless of how big/small of a time step
% that we choose.

%% Problem 3(c)

clear; close all; clc;

% Parameters
dt1 = 1;
dt2 = 0.5;
dt3 = 0.25;
T = 10;
func = @odefun_3c;

% Solve the ivp, odefun
[tt1,uu1] = crank_nicolson(func,T,dt1);
[tt2,uu2] = crank_nicolson(func,T,dt2);
[tt3,uu3] = crank_nicolson(func,T,dt3);

% Plot the solutions
figure(3)
plot(tt1,uu1);
hold on
plot(tt2,uu2);
plot(tt3,uu3);
syms t
ezplot(cos(t),[0,10]);
xlabel('time'); 
ylabel('u');
title('Q3(c): Crank-Nicolson solutions');
legend('dt=1','dt=0.5','dt=0.25','analytical')
hold off

%% Problem 3(e)

% Consider the IVP (1) introduced in Problem 2 with ? = 1000. 
% Compute the error at the ?nal time |u(T)? ˜ u(T)| for the Crank-Nicolson
% method with time steps of 1, 0.5, 0.25, 0.125, 0.0625. Does the observed 
% convergence rate match your expectation?

clear; close all; clc;

% Parameters
dt1 = 1;
dt2 = 0.5;
dt3 = 0.25;
dt4 = 0.125;
dt5 = 0.0625;
T = 10;
func = @odefun_3c;

% Solve the ivp, odefun
[tt1,uu1] = crank_nicolson(func,T,dt1);
[tt2,uu2] = crank_nicolson(func,T,dt2);
[tt3,uu3] = crank_nicolson(func,T,dt3);
[tt4,uu4] = crank_nicolson(func,T,dt4);
[tt5,uu5] = crank_nicolson(func,T,dt5);

% Compute the errors
err1 = abs(cos(10) - uu1(end));
err2 = abs(cos(10) - uu2(end));
err3 = abs(cos(10) - uu3(end));
err4 = abs(cos(10) - uu4(end));
err5 = abs(cos(10) - uu5(end));

% log-log plot of error vs. timestep
dts = [1, 0.5, 0.25, 0.125,0.0625]';
errs = [err1, err2, err3, err4, err5]';
figure(4)
loglog(dts,errs,'-o','LineWidth',1.5);
title('Q3(e): Convergence rate of Crank Nicolson');
xlabel('dt');
ylabel('error at T=10');
hold off

% convergence rate of FE
p = polyfit(log(dts),log(errs),1);
conv_rate_CN = p(1)

% COMMENTS: We observe a convergence rate of approximately 2 (precisely
% 2.0329). This is expected based on our previous analysis, which shows
% that the order of accuracy of the Crank Nicolson scheme is in fact 2.


%% Problem 2(b)

function [f,dfdu] = odefun(t,u)
    gamma = 1;
    f = gamma * (cos(t) - u) - sin(t);
    dfdu = -gamma;
end

%% Problem 2(c)

function [tt,uu] = forward_euler(odefun,T,dt)
    uu(1) = 1; % initial condition
    tt = (0:dt:T)';
    N = length(tt);
    for i = 2:N
        [f,dfdu] = odefun(tt(i-1),uu(i-1));
        uu(i) = uu(i-1) + dt*f;
    end
end

%% Problem 3(a)

function [tt,uu] = crank_nicolson(odefun,T,dt)
    tol = 10^-12; % we set the tol arbitrarily
    uu(1) = 1; % initial condition
    tt = (0:dt:T)';
    N = length(tt);
    for j = 2:N
        % perform Newton's method at each time step for uj
        uu(j) = uu(j-1); % let the initial guess be previous time step
        [f_prev,dfdu_prev] = odefun(tt(j-1),uu(j-1));
        [f_curr,dfdu_curr] = odefun(tt(j),uu(j));
        err = abs(-0.5*dt*(f_prev + f_curr));
        while err > tol
            r = uu(j) - uu(j-1) - 0.5*dt*f_curr - 0.5*dt*f_prev;
            dr = 1 - 0.5*dt*dfdu_curr;
            uu(j) = uu(j) - r/dr;
            [f_curr,dfdu_curr] = odefun(tt(j),uu(j));
            err = abs(uu(j) - uu(j-1) - 0.5*dt*(f_prev + f_curr));
        end
    end
end


%% Other function definitinos 

% ode function for Problem 3(b)

function [f,dfdu] = odefun_3b(t,u)
    f = 2*(t^3 + t)/u;
    dfdu = -2*(t^3 + t)/u^2;
end

% ode function for Problem 3(c): gamma=1000

function [f,dfdu] = odefun_3c(t,u)
    gamma = 1000;
    f = gamma * (cos(t) - u) - sin(t);
    dfdu = -gamma;
end
