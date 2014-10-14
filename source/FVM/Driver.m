%% Driver: produce the results for the SEB410 assignment.
% This script generates all the results included in Morgan Wall's report for
% the SEB410 assignment. The script contains solutions for three (3) problems
% aimed to test a vertex-centred finite volume method for solving generalised
% two-dimensional non-linear advection-diffusion equations.

clear all;
close all;

%% Problem 1: Linear convection-diffusion (with analytic solution)

% Initialise problem parameters
tFinal = 1.25;
Dxx = @(phi) ones(length(phi), 1) * 0.01;
Dyy = @(phi) ones(length(phi), 1) * 0.01;
Vx = @(phi) ones(length(phi), 1) .* 0.8;
Vy = @(phi) ones(length(phi), 1) .* 0.8;
source = @(phi) phi .* 0;

% Initialise boundary conditions
northBC = struct('A', 0, 'B', 1, 'C', 0);
eastBC = struct('A', 0, 'B', 1, 'C', 0);
southBC = struct('A', 0, 'B', 1, 'C', 0);
westBC = struct('A', 0, 'B', 1, 'C', 0);

% Construct initial condition
initialCondition = zeros(length(nodesY), length(nodesX));
initialCondition(round(length(nodesY) / 2), round(length(nodesY) / 2)) = 1;

% Construct mesh
xLower = 0;
xUpper = 2;
xCount = 30;
xGeoParameters = struct('lowerIsGeometric', false, ...
    'upperIsGeometric', false, 'commonRatio', 1); 

yLower = 0;
yUpper = 2;
yCount = 30;
yGeoParameters = struct('lowerIsGeometric', false, ...
    'upperIsGeometric', false, 'commonRatio', 1);

[nodesX, nodesY] = GenerateStructuredMesh(xLower, xUpper, xCount, ...
    yLower, yUpper, yCount, xGeoParameters, yGeoParameters);
nodesY = flipud(nodesY);

rows = length(nodesY);
columns = length(nodesX);

% Initialise solver parameters
theta = 1;
advectionHandling = 'averaging';

dt = 0.001;
storedTimeSteps = 100;


newtonParameters = struct('maxIterations', 5, 'tolUpdate', 1e-6, ...
    'tolResidual', 1e-6);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 20, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'none', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

% Solve problem
[tout, yout] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters);

% Output plots and metrics
figure;

surf(nodesX, nodesY, reshape(yout(:, end), rows, columns));
plotTitle = ['Test Problem (G1.1): Gaussian Diffusion (t = ' ...
    num2str(tout(end)) ')'];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

figure;

surf(nodesX, nodesY, reshape(yout(:, 1), rows, columns));
plotTitle = 'Test Problem (G1.2): Gaussian Diffusion (t = 0)';
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

%% Problem 2: Non-linear convection-diffusion (with known steady-state)

%% Problem 3: Non-linear convection-diffusion (incl. Peclet number)


