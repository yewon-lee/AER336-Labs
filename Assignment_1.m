%% AER336 Assigment 1

% Yewon Lee
% Student number: 1004407260

%% Problem 2 (b)

M = 1.3;
p = mach2pressure_pw1(M)

%% Problem 2 (c)

% original data set
M_orig = linspace(0,4,9);
p_orig = [
1.000000000000000
0.843019175422553
0.528281787717174
0.272403066476657
0.127804525462951
0.058527663465935
0.027223683703863
0.013110919994476
0.006586087307989]';

% evaluate the interpolant
M_interp = linspace(0,4,8001);
p_interp = NaN(1,length(M_interp));
for i = 1:length(M_interp)-1 
    p_interp(i) = mach2pressure_pw1(M_interp(i));
end
p_interp(length(p_interp)) = p_orig(length(p_orig));

% plot interpolant and data points
figure(1)
plot(M_interp,p_interp,'LineWidth',2);
hold on
plot(M_orig,p_orig,'o','LineWidth',2);
xlabel('M');
ylabel('p/p_0');
legend('Interpolant','Given data');
hold off

%% Problem 2 (d)

% analytical form of the pressure ratio vs. Mach number
gamma = 1.4;
p_true = @(M) (1 + (gamma-1)/2 * M^2)^(-gamma/(gamma - 1));

% evaluate the interpolant
M_interp = linspace(0,4,8001);
p_interp = NaN(1,length(M_interp));
for i = 1:length(M_interp)-1 
    p_interp(i) = mach2pressure_pw1(M_interp(i));
end
p_interp(length(p_interp)) = 0.006586087307989;

% evaluate p_true_vec
p_true_vec = NaN(1,length(M_interp));
for i = 1:length(p_true_vec)
    p_true_vec(i) = p_true(M_interp(i));
end

% error vector
err = abs(p_interp - p_true_vec);

% plot error vector as a function of M
figure(2)
plot(M_interp,err,'LineWidth',2);
xlabel('M');
ylabel('error');
hold off

%% Problem 3 (b)

% Here we construct and output the re-ordered x and z vectors

% Load the original data
load('rae2822geom.dat');
x_orig = rae2822geom(:,1);
zl_orig = rae2822geom(:,2);
zu_orig = rae2822geom(:,3);

% Construct re-ordered x and z vectors
x = NaN(33,1); z = NaN(33,1);
x(1:17) = x_orig(17:-1:1);
x(18:33) = x_orig(2:17);
z(1:17) = zu_orig(17:-1:1);
z(18:33) = zl_orig(2:17);

% Print the re-ordered x and z vectors
x
z

%% Problem 3 (c): Part 2 (evaluating li for the last interpolation point)

% Note: for the function definition, please refer to Part 1 of Problem 3
% (c), defined later in this file

% Load the original data
load('rae2822geom.dat');
x_orig = rae2822geom(:,1);
zl_orig = rae2822geom(:,2);
zu_orig = rae2822geom(:,3);

% Construct re-ordered x and z vectors
x = NaN(33,1); z = NaN(33,1);
x(1:17) = x_orig(17:-1:1);
x(18:33) = x_orig(2:17);
z(1:17) = zu_orig(17:-1:1);
z(18:33) = zl_orig(2:17);

% Compute the li vector as well as li(33)
li = arclength(x,z);
li(33)

%% Problem 3 (d)

% Load the original data
load('rae2822geom.dat');
x_orig = rae2822geom(:,1);
zl_orig = rae2822geom(:,2);
zu_orig = rae2822geom(:,3);

% Construct re-ordered x and z vectors
x = NaN(33,1); z = NaN(33,1);
x(1:17) = x_orig(17:-1:1);
x(18:33) = x_orig(2:17);
z(1:17) = zu_orig(17:-1:1);
z(18:33) = zl_orig(2:17);

% Compute the li vector
li = arclength(x,z);

% Use spline to find x_spline and y_spline as functions of l
l_spline = NaN(7999,1);
l_spline(1:4000) = linspace(li(1),li(17),4000)';
diff = abs(l_spline(1) - l_spline(2));
l_spline(4001:7999) = linspace(li(17) + diff, li(end), 3999)';
x_spline = spline(li, x, l_spline);
y_spline = spline(li, z, l_spline);

% Plot the cubic spline
figure(4)
plot(x_spline, y_spline,'LineWidth',2);
hold on
plot(x,z,'o','LineWidth',1);
legend('cubic spline approximation','given data')
axis([0 1 -0.5 0.5])
title('RAE2822 airfoil')
hold off

% focus on the leading edge
figure(5)
plot(x_spline, y_spline,'LineWidth',2);
hold on
plot(x,z,'o','LineWidth',1);
legend('cubic spline approximation','given data')
axis([-0.02 0.02 -0.02 0.02])
title('RAE2822 airfoil')
hold off

% COMMENT: We can see from figure 5 (the figure that focuses on the leading 
% edge of the airfoil) that the cubic spline approximation is indeed smooth
% at the leading edge

%% Problem 2 (a)

% Piecewise linear interpolation function

function p = mach2pressure_pw1(M)

% initialize M_vec and p_vec
M_vec = linspace(0,4,9);
p_vec = [
1.000000000000000
0.843019175422553
0.528281787717174
0.272403066476657
0.127804525462951
0.058527663465935
0.027223683703863
0.013110919994476
0.006586087307989];

% locate segment to which a given M belongs
k = (M - mod(M,0.5))/0.5 + 1;

% perform linear interpolation in segment Sk
p = p_vec(k) + ((p_vec(k+1) - p_vec(k)) / 0.5)*(M - M_vec(k));

end

%% Problem 3 (c): Part 1 (function definition)

function li = arclength(xi,yi)

li = NaN(length(xi),1);
for k = 1:length(xi)
    if k == 1
        li(k) = 0;
    else
        li(k) = li(k-1) + sqrt( (xi(k) - xi(k-1))^2 + (yi(k) - yi(k-1))^2 );
    end 
end

end


