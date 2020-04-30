clear;
close all;
clc;

%% Initialization
COLLIMATOR_SIDE = 4;
EX = 10; % Ellipse major axis range
EY = 5; % Ellipse minor axis range
EZ = 2; % Ellipse z axis range

%% Create symbols and ellipse equations
syms x y z E T O1 O2 O3;
E = x^2/(EX^2) + y^2/(EY^2) - 1;
T = x^2 + y^2 + z^2 - 1;
O1 = (x-3)^2 + 4*y^2 + z^2/4 - 1;
O2 = 4*(x-3)^2 + 4*(y-3)^2 + 4*z^2 - 1;
O3 = 4*x^2 + 9*(y-2)^2 + 16*z^2 - 1;

%% Plot tumor, critical organs, normal tissue
figure;
hold on;
grid on;
ezplot(E, [-1*EX, EX, -1 * EY, EY])
ezplot(subs(T, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O1, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O2, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O3, z, 0), [-1*EX, EX, -1 * EY, EY])
xlim([-1*(EX+0.5), EX+0.5])
ylim([-1*(EY+0.5), EY+0.5])

%% Create the mesh
% Make the precision in the ellipse is smaller than that of the collimator.
N = 20;
Nx = N + 1;
Ny = N + 1;
Nz = N + 1;
dx = 2*EX/Nx; % Range along x-axis is 20
dy = 2*EY/Ny; % Range along y-axis is 10
dz = 2*EZ/Nz; % Range along z-axis is 4

x_pts = linspace(-1*EX, EX, Nx);
y_pts = linspace(-1*EY, EY, Nx);
z_pts = linspace(-1*EZ, EZ, Nx);
[X, Y, Z] = meshgrid(x_pts, y_pts, z_pts);
X = X(:); Y = Y(:); Z = Z(:);
P = [X, Y, Z];
% Plot the mesh
scatter(X, Y, 'b')

% Find the points that fall inside the ellipse
in_eps_flags = double(subs(E, {x, y}, {P(:,1), P(:,2)})) <= 0;

% Plot the points that fall inside the ellipse
Q = P(in_eps_flags, :);
scatter(Q(:, 1), Q(:, 2), 'r')

%% Discretize the collimator.
% It is in the Y-Z plane
zc_pts = linspace(-1*EZ, EZ, 8+1); 
yc_pts = linspace(-1*EZ, EZ, 8+1); 
[Z, Y] = meshgrid(zc_pts, yc_pts);
Z = Z(:); Y = Y(:);

%% Plot tumor, critical organs, normal tissue
E = x^2/(EX^2) + y^2/(EY^2) - 1;
T = x^2 + y^2 + z^2 - 1;
O1 = (x-3)^2 + 4*y^2 + z^2/4 - 1;
O2 = 4*(x-3)^2 + 4*(y-3)^2 + 4*z^2 - 1;
O3 = 4*x^2 + 9*(y-2)^2 + 16*z^2 - 1;
O4 = (x+3)^2 + 4*y^2 + z^2/4 - 1;
O5 = 4*(x+3)^2 + 4*(y-3)^2 + 4*z^2 - 1;
O6 = 4*x^2 + 9*(y+2)^2 + 16*z^2 - 1;
O7 = 4*(x+3)^2 + 4*(y+3)^2 + 4*z^2 - 1;
O8 = 4*(x-3)^2 + 4*(y+3)^2 + 4*z^2 - 1;

figure;
hold on;
grid on;
ezplot(E, [-1*EX, EX, -1 * EY, EY])
ezplot(subs(T, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O1, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O2, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O3, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O4, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O5, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O6, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O7, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O8, z, 0), [-1*EX, EX, -1 * EY, EY])
xlim([-1*(EX+0.5), EX+0.5])
ylim([-1*(EY+0.5), EY+0.5])
xlabel('X-Axis', 'FontSize', 14)
ylabel('Y-Axis', 'FontSize', 14)
title('', 'FontSize', 18)
