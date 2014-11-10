%% Driver: produce the results for the SEB410 assignment.
% This script generates all the results included in Morgan Wall's report for
% the SEB410 assignment. The script contains solutions for three (3) problems
% aimed to test a vertex-centred finite volume method for solving generalised
% two-dimensional non-linear advection-diffusion equations.

clear all;
close all;

%% Problem 1: Linear convection-diffusion (with analytic solution)

% Initialise problem parameters
dt = 0.001;
tFinal = 1.25;

xC = 0.5;
yC = 0.5;

DXX = 0.01;
Dxx = @(phi) phi .* 0 + DXX;

DYY = 0.01;
Dyy = @(phi) phi .* 0 + DYY;

VX = 0.8;
Vx = @(phi) phi .* 0 + VX;

VY = 0.8;
Vy = @(phi) phi .* 0 + VY;

source = @(phi) phi .* 0;

% Construct mesh
xLower = 0;
xUpper = 2;
xCount = 85;
xGeoParameters = struct('lowerIsGeometric', false, ...
    'upperIsGeometric', false, 'commonRatio', 1); 

yLower = 0;
yUpper = 2;
yCount = 85;
yGeoParameters = struct('lowerIsGeometric', false, ...
    'upperIsGeometric', false, 'commonRatio', 1);

[nodesX, nodesY] = GenerateStructuredMesh(xLower, xUpper, xCount, ...
    yLower, yUpper, yCount, xGeoParameters, yGeoParameters);
nodesY = flipud(nodesY);

rows = length(nodesY);
columns = length(nodesX);

% Initialise analytic solution
phiAnalytic = @(x, y, t) exp( -(x - VX * t - xC).^2 ./ (DXX * (4 * t + 1)) ...
    - (y - VY * t - yC).^2 ./ (DYY * (4 * t + 1)) ) ./(4 * t + 1);

% Initialise boundary conditions
dirichletHackCoef = 10000;

northC = @(x, t) dirichletHackCoef .* phiAnalytic(x, yUpper, t);
northBC = struct('A', dirichletHackCoef, 'B', 1, 'C', northC);

eastC = @(y, t) dirichletHackCoef .* phiAnalytic(xUpper, y, t);
eastBC = struct('A', dirichletHackCoef, 'B', 1, 'C', eastC);

southC = @(x, t) dirichletHackCoef .* phiAnalytic(x, yLower, t);
southBC = struct('A', dirichletHackCoef, 'B', 1, 'C', southC);

westC = @(y, t) dirichletHackCoef .* phiAnalytic(xLower, y, t);
westBC = struct('A', dirichletHackCoef, 'B', 1, 'C', westC);

% Construct initial condition
[X, Y] = meshgrid(nodesX(:), nodesY(:));
initialCondition = phiAnalytic(X(:), Y(:), 0);

% Construct analytic solution (at final time)
analyticSolution = phiAnalytic(X(:), Y(:), tFinal);

% Initialise solver parameters
theta = 1;
advectionHandling = 'upwinding';

storedTimeSteps = 100;

newtonParameters = struct('maxIterations', 5, 'tolUpdate', 1e-8, ...
    'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-6, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'none', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

% Solve problem
tic;
[tout, yout] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters);
toc;

% Output plots and metrics
figure;

surf(nodesX, nodesY, reshape(yout(:, 1), rows, columns));
plotTitle = 'Analytic Solution (Problem 1) where t = 0';
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

figure;

subplot(2, 1, 1);
surf(nodesX, nodesY, reshape(yout(:, end), rows, columns));
plotTitle = ['Numeric Solution (Problem 1) where t = ' num2str(tFinal)];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

subplot(2, 1, 2);
surf(nodesX, nodesY, reshape(analyticSolution(:, end), rows, columns));
plotTitle = ['Analytic Solution (Problem 1) where t = ' num2str(tFinal)];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

error = norm(yout(:, end) - analyticSolution(:, end)) / length(yout(:, end));

%% Problem 2: Non-linear convection-diffusion (with known steady-state)

%% Problem 3: Non-linear convection-diffusion (incl. Peclet number)
