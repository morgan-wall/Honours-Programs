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

source = @(x, y) x .* 0;

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

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

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

surf(nodesX, nodesY, reshape(yout(:, end), rows, columns));
plotTitle = ['Numeric Solution (Problem 1) where t = ' num2str(tFinal)];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

figure;

surf(nodesX, nodesY, reshape(analyticSolution(:, end), rows, columns));
plotTitle = ['Analytic Solution (Problem 1) where t = ' num2str(tFinal)];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

error = norm(yout(:, end) - analyticSolution(:, end)) / sqrt(length(yout(:, end)));

%% Problem 2: Non-linear convection-diffusion (with known steady-state)

% Initialise problem parameters
dt = 0.001;
tFinal = 10;

Dxx = @(phi) phi .^ 2;

Dyy = @(phi) 5 .* phi .^ 2;

C = 20;
Vx = @(phi) (C / 2) .* phi;
Vy = @(phi) (C / 2) .* phi;

source = @(x,y) 100 * x .* y .* (1-x) .* (1-y) .* exp(2 * x.^4.5) ...
    .* (C * y .* (1-y) .* (1 - 2 * x + 4.5 * (1-x) .* x.^4.5) ...
    + C * x .* (1-x) .* (1 - 2 * y) - 20 * y.^2 .* (1-y).^2 .* exp(x.^4.5) ...
    .* (1 + 4 * x .^ 2 + 4.5^2 * x.^9 .* (1 + x.^2 - 2 * x) - 4 * x ...
    + 9 * x.^4.5 .* (1 - 3 * x + 2 * x.^2)) - 100 * exp(x.^4.5) .* x.^2 ...
    .* (1-x).^2 .* (1 - 2 * y).^2 - 10 * x .* y.^2 .* (1-x) .* (1-y).^2 ...
    .* exp(x.^4.5) .* (4.5 * x.^3.5 .*(5.5 - 7.5 * x) + 4.5^2 * x.^8 ...
    .* (1-x) - 2) + 100 * x.^ 2 .* y .* (1-x).^2 .* (1-y) .* exp(x.^4.5));

% sourceMaple = @(x, y) 0.2000e4 .* x .* (y .^ 2) .* ((1 - x) .^ 2) .* ((1 - y) .^ 2) ...
%     .* exp((x .^ 0.45e1)) .^ 2 - 0.2000e4 .* (x .^ 2) .* (y .^ 2) .* (1 - x) ...
%     .* ((1 - y) .^ 2) .* exp((x .^ 0.45e1)) .^ 2 + 0.90000e4 .* (x .^ 0.55e1) ...
%     .* (y .^ 2) .* ((1 - x) .^ 2) .* ((1 - y) .^ 2) .* exp((x .^ 0.45e1)) .^ 2 ...
%     - 0.200e3 .* x .* (y .^ 2) .* ((1 - x) .^ 2) .* ((1 - y) .^ 2) ...
%     .* exp((x .^ 0.45e1)) .^ 2 .* (0.10e2 .* y .* (1 - x) .* (1 - y) ...
%     .* exp((x .^ 0.45e1)) - 0.10e2 .* x .* y .* (1 - y) .* exp((x .^ 0.45e1)) ...
%     + 0.450e2 .* (x .^ 0.45e1) .* y .* (1 - x) .* (1 - y) .* exp((x .^ 0.45e1))) ...
%     + 0.200e3 .* (x .^ 2) .* (y .^ 2) .* (1 - x) .* ((1 - y) .^ 2) ...
%     .* exp((x .^ 0.45e1)) .^ 2 .* (0.10e2 .* y .* (1 - x) .* (1 - y) ...
%     .* exp((x .^ 0.45e1)) - 0.10e2 .* x .* y .* (1 - y) .* exp((x .^ 0.45e1)) ...
%     + 0.450e2 .* (x .^ 0.45e1) .* y .* (1 - x) .* (1 - y) .* exp((x .^ 0.45e1))) ...
%     - 0.9000e3 .* (x .^ 0.55e1) .* (y .^ 2) .* ((1 - x) .^ 2) .* ((1 - y) .^ 2) ...
%     .* exp((x .^ 0.45e1)) .^ 2 .* (0.10e2 .* y .* (1 - x) .* (1 - y) ...
%     .* exp((x .^ 0.45e1)) - 0.10e2 .* x .* y .* (1 - y) .* exp((x .^ 0.45e1)) ...
%     + 0.450e2 .* (x .^ 0.45e1) .* y .* (1 - x) .* (1 - y) .* exp((x .^ 0.45e1))) ...
%     - 0.100e3 .* (x .^ 2) .* (y .^ 2) .* ((1 - x) .^ 2) .* ((1 - y) .^ 2) ...
%     .* exp((x .^ 0.45e1)) .^ 2 .* (-0.20e2 .* y .* (1 - y) .* exp((x .^ 0.45e1)) ...
%     + 0.24750e3 .* y .* (1 - x) .* (1 - y) .* (x .^ 0.35e1) .* exp((x .^ 0.45e1)) ...
%     - 0.900e2 .* (x .^ 0.45e1) .* y .* (1 - y) .* exp((x .^ 0.45e1)) + 0.20250e3 ...
%     .* (x .^ 0.80e1) .* y .* (1 - x) .* (1 - y) .* exp((x .^ 0.45e1))) + 0.2000e4 ...
%     .* (x .^ 2) .* y .* ((1 - x) .^ 2) .* ((1 - y) .^ 2) .* exp((x .^ 0.45e1)) .^ 2 ...
%     - 0.2000e4 .* (x .^ 2) .* (y .^ 2) .* ((1 - x) .^ 2) .* (1 - y) ...
%     .* exp((x .^ 0.45e1)) .^ 2 - 0.1000e4 .* (x .^ 2) .* y .* ((1 - x) .^ 2) ...
%     .* ((1 - y) .^ 2) .* exp((x .^ 0.45e1)) .^ 2 .* (0.10e2 .* x .* (1 - x) .* (1 - y) ...
%     .* exp((x .^ 0.45e1)) - 0.10e2 .* x .* y .* (1 - x) .* exp((x .^ 0.45e1))) ...
%     + 0.1000e4 .* (x .^ 2) .* (y .^ 2) .* ((1 - x) .^ 2) .* (1 - y) ...
%     .* exp((x .^ 0.45e1)) .^ 2 .* (0.10e2 .* x .* (1 - x) .* (1 - y) ...
%     .* exp((x .^ 0.45e1)) - 0.10e2 .* x .* y .* (1 - x) .* exp((x .^ 0.45e1))) ...
%     + 0.10000e5 .* (x .^ 3) .* (y .^ 2) .* ((1 - x) .^ 3) .* ((1 - y) .^ 2) ...
%     .* exp((x .^ 0.45e1)) .^ 3;

% Construct mesh
xLower = 0;
xUpper = 1;
xCount = 100;
xGeoParameters = struct('lowerIsGeometric', false, ...
    'upperIsGeometric', false, 'commonRatio', 1); 

yLower = 0;
yUpper = 1;
yCount = 100;
yGeoParameters = struct('lowerIsGeometric', false, ...
    'upperIsGeometric', false, 'commonRatio', 1);

[nodesX, nodesY] = GenerateStructuredMesh(xLower, xUpper, xCount, ...
    yLower, yUpper, yCount, xGeoParameters, yGeoParameters);
nodesY = flipud(nodesY);

rows = length(nodesY);
columns = length(nodesX);

% Determine source term
[X, Y] = meshgrid(nodesX(:), nodesY(:));
sourceTerm = source(X(:), Y(:));
% sourceTermMaple = sourceMaple(X(:), Y(:));

% Initialise boundary conditions
dirichletHackCoef = 10000;

northC = @(x, t) zeros(length(x), 1);
northBC = struct('A', dirichletHackCoef, 'B', 1, 'C', northC);

eastC = @(y, t) zeros(length(y), 1);
eastBC = struct('A', dirichletHackCoef, 'B', 1, 'C', eastC);

southC = @(x, t) zeros(length(x), 1);
southBC = struct('A', dirichletHackCoef, 'B', 1, 'C', southC);

westC = @(y, t) zeros(length(y), 1);
westBC = struct('A', dirichletHackCoef, 'B', 1, 'C', westC);

% Construct initial condition
initialCondition = zeros(rows, columns);

% Construct steady state solution
steadyState = @(x, y) 10 * x .* y .* (1 - x) .* (1 - y) .* exp(x .^ 4.5);
steadyStateSolution = steadyState(X(:), Y(:));

% Initialise solver parameters
theta = 1;
advectionHandling = 'averaging';

storedTimeSteps = 1000;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 5000, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

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
plotTitle = 'Analytic Solution (Problem 2) where t = 0';
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

figure;
surf(nodesX, nodesY,reshape(sourceTerm, rows, columns));
plotTitle = 'Source Term (Problem 2)';
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

figure;

surf(nodesX, nodesY, reshape(yout(:, end), rows, columns));
plotTitle = ['Numeric Solution (Problem 2) where t = ' num2str(tFinal)];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

figure;

surf(nodesX, nodesY, reshape(steadyStateSolution, rows, columns));
plotTitle = 'Steady State Solution (Problem 2)';
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

%% Problem 3: Non-linear convection-diffusion (incl. Peclet number)
