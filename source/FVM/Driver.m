%% Driver: unit tests for the Solver.m file.
% This script is used for testing purposes. The script contains a series of
% "unit" tests for ensuring the the finite-volume solver referenced in 
% Solver.m is behaving correctly. Strictly speaking, this script does not
% contain unit tests, rather it conducts a series of simulations of various
% non-linear, two-dimensional advection-diffusion equations. This script
% aims to promote test driven development.
%

clear all;
close all;

% Initialise solver parameters
theta = 0;
advectionHandling = 'averaging';

%% Test: Gaussian Diffusion (G1) - Single Point Initial Condition

% Initialise temporal parameters
tFinal = 0.1;
storedTimeSteps = 100;
% tFinal = 0.0005;
% storedTimeSteps = 1;

% Initialise equation parameters
Dxx = @(phi) 0.1;
Dyy = @(phi) 0.1;
Vx = @(phi) 0;
Vy = @(phi) 0;
source = @(phi) 0;

% Initialise mesh parameters
nodesX = 0:0.05:1;
nodesY = 1:-0.05:0;

rows = length(nodesY);
columns = length(nodesX);

% Initialise boundary conditions
northBC = struct('A', 0, 'B', 1, 'C', 0);
eastBC = struct('A', 0, 'B', 1, 'C', 0);
southBC = struct('A', 0, 'B', 1, 'C', 0);
westBC = struct('A', 0, 'B', 1, 'C', 0);

% Construct initial condition
initialCondition = zeros(length(nodesY), length(nodesX));
initialCondition(round(length(nodesY) / 2), round(length(nodesY) / 2)) = 5;

% Solve problem
[tout, yout] = Solver(tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps);

% Output plots and metrics
figure;

surf(nodesX, nodesY, reshape(yout(:, end), rows, columns));
plotTitle = ['Test Problem (G1.1): Gaussian Diffusion (t = ' ...
    num2str(tout(end)) ')'];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

disp(['Total at end time: ' num2str(sum(sum(yout(:, end)))) '.']);

figure;

surf(nodesX, nodesY, reshape(yout(:, 1), rows, columns));
plotTitle = 'Test Problem (G1.2): Gaussian Diffusion (t = 0)';
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

%% Test: Dirichlet & Neumann Boundary Conditions (N1) - Input at West face. 

% Initialise temporal parameters
tFinal = 0.1;
storedTimeSteps = 100;

% Initialise equation parameters
Dxx = @(phi) 0.1;
Dyy = @(phi) 0.1;
Vx = @(phi) 0;
Vy = @(phi) 0;
source = @(phi) 0;

% Initialise mesh parameters
nodesX = 0:0.05:1;
nodesY = 1:-0.05:0;

rows = length(nodesY);
columns = length(nodesX);

% Initialise boundary conditions
northBC = struct('A', 0, 'B', 1, 'C', 0);
eastBC = struct('A', 0, 'B', 1, 'C', 0);
southBC = struct('A', 0, 'B', 1, 'C', 0);
westBC = struct('A', 1000, 'B', 1, 'C', 1000);

% Construct initial condition
initialCondition = zeros(length(nodesY), length(nodesX));
initialCondition(:, 1) = 1;

% Solve problem
[tout, yout] = Solver(tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps);

% Output plots and metrics
figure;

surf(nodesX, nodesY, reshape(yout(:, end), rows, columns));
plotTitle = ['Test Problem (N1.1): Dirichlet & Neumann Boundary ' ...
    'Conditions (t = ' num2str(tout(end)) ')'];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

figure;

surf(nodesX, nodesY, reshape(yout(:, 1), rows, columns));
plotTitle = ['Test Problem (N1.2): Dirichlet & Neumann Boundary ' ...
    'Conditions (t = 0 )'];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

%% Test: Dirichlet & Neumann Boundary Conditions (N2) - Input at North face. 

% Initialise temporal parameters
tFinal = 0.1;
storedTimeSteps = 100;

% Initialise equation parameters
Dxx = @(phi) 0.1;
Dyy = @(phi) 0.1;
Vx = @(phi) 0;
Vy = @(phi) 0;
source = @(phi) 0;

% Initialise mesh parameters
nodesX = 0:0.05:1;
nodesY = 1:-0.05:0;

rows = length(nodesY);
columns = length(nodesX);

% Initialise boundary conditions
northBC = struct('A', 1000, 'B', 1, 'C', 1000);
eastBC = struct('A', 0, 'B', 1, 'C', 0);
southBC = struct('A', 0, 'B', 1, 'C', 0);
westBC = struct('A', 0, 'B', 1, 'C', 0);

% Construct initial condition
initialCondition = zeros(length(nodesY), length(nodesX));
initialCondition(1, :) = 1;

% Solve problem
[tout, yout] = Solver(tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps);

% Output plots and metrics
figure;

surf(nodesX, nodesY, reshape(yout(:, end), rows, columns));
plotTitle = ['Test Problem (N2.1): Dirichlet & Neumann Boundary ' ...
    'Conditions (t = ' num2str(tout(end)) ')'];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

figure;

surf(nodesX, nodesY, reshape(yout(:, 1), rows, columns));
plotTitle = ['Test Problem (N2.2): Dirichlet & Neumann Boundary ' ...
    'Conditions (t = 0 )'];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

%% Test: Dirichlet & Neumann Boundary Conditions (N3) - Input at East face. 

% Initialise temporal parameters
tFinal = 0.1;
storedTimeSteps = 100;

% Initialise equation parameters
Dxx = @(phi) 0.1;
Dyy = @(phi) 0.1;
Vx = @(phi) 0;
Vy = @(phi) 0;
source = @(phi) 0;

% Initialise mesh parameters
nodesX = 0:0.05:1;
nodesY = 1:-0.05:0;

rows = length(nodesY);
columns = length(nodesX);

% Initialise boundary conditions
northBC = struct('A', 0, 'B', 1, 'C', 0);
eastBC = struct('A', 1000, 'B', 1, 'C', 1000);
southBC = struct('A', 0, 'B', 1, 'C', 0);
westBC = struct('A', 0, 'B', 1, 'C', 0);

% Construct initial condition
initialCondition = zeros(length(nodesY), length(nodesX));
initialCondition(:, end) = 1;

% Solve problem
[tout, yout] = Solver(tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps);

% Output plots and metrics
figure;

surf(nodesX, nodesY, reshape(yout(:, end), rows, columns));
plotTitle = ['Test Problem (N3.1): Dirichlet & Neumann Boundary ' ...
    'Conditions (t = ' num2str(tout(end)) ')'];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

figure;

surf(nodesX, nodesY, reshape(yout(:, 1), rows, columns));
plotTitle = ['Test Problem (N3.2): Dirichlet & Neumann Boundary ' ...
    'Conditions (t = 0 )'];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

%% Test: Dirichlet & Neumann Boundary Conditions (N4) - Input at South face. 

% Initialise temporal parameters
tFinal = 0.1;
storedTimeSteps = 100;

% Initialise equation parameters
Dxx = @(phi) 0.1;
Dyy = @(phi) 0.1;
Vx = @(phi) 0;
Vy = @(phi) 0;
source = @(phi) 0;

% Initialise mesh parameters
nodesX = 0:0.05:1;
nodesY = 1:-0.05:0;

rows = length(nodesY);
columns = length(nodesX);

% Initialise boundary conditions
northBC = struct('A', 0, 'B', 1, 'C', 0);
eastBC = struct('A', 0, 'B', 1, 'C', 0);
southBC = struct('A', 1000, 'B', 1, 'C', 1000);
westBC = struct('A', 0, 'B', 1, 'C', 0);

% Construct initial condition
initialCondition = zeros(length(nodesY), length(nodesX));
initialCondition(end, :) = 1;

% Solve problem
[tout, yout] = Solver(tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps);

% Output plots and metrics
figure;

surf(nodesX, nodesY, reshape(yout(:, end), rows, columns));
plotTitle = ['Test Problem (N4.1): Dirichlet & Neumann Boundary '...
    'Conditions (t = ' num2str(tout(end)) ')'];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

figure;

surf(nodesX, nodesY, reshape(yout(:, 1), rows, columns));
plotTitle = ['Test Problem (N4.2): Dirichlet & Neumann Boundary ' ...
    'Conditions (t = 0 )'];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');
