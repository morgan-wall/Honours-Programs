%% Driver_p3: produce the results for the SEB410 assignment.
% This script generates all the results included in Morgan Wall's report for
% the SEB410 assignment. The script contains solutions for the third problem
% aimed to test a vertex-centred finite volume method for solving generalised
% two-dimensional non-linear advection-diffusion equations.

clear all;
close all;

% Initialise problem parameters
dt = 0.0001;
tFinal = 7;

Pe = 10;

DXX = 1 / Pe;
Dxx = @(phi, x, y, t) phi .* 0 + DXX;

DYY = 1 / Pe;
Dyy = @(phi, x, y, t) phi .* 0 + DYY;

Vx = @(phi, x, y, t) 2 .* y .* (1 - x.^2);

Vy = @(phi, x, y, t) -2 .* x .* (1 - y.^2);

source = @(phi, x, y, t) x .* 0;

% Construct mesh
xLower = -1;
xUpper = 1;
xCount = 201;
xGeoParameters = struct('lowerIsGeometric', false, ...
    'upperIsGeometric', false, 'commonRatio', 1); 

yLower = 0;
yUpper = 1;
yCount = 101;
yGeoParameters = struct('lowerIsGeometric', false, ...
    'upperIsGeometric', false, 'commonRatio', 1);

[nodesX, nodesY] = GenerateStructuredMesh(xLower, xUpper, xCount, ...
    yLower, yUpper, yCount, xGeoParameters, yGeoParameters);
nodesY = flipud(nodesY);

rows = length(nodesY);
columns = length(nodesX);

% Initialise boundary conditions
dirichletHackCoef = 10000;

northA = @(x, t) x .* 0 + dirichletHackCoef;
northB = @(x, t) x .* 0 + 1;
northC = @(x, t) dirichletHackCoef .* ones(length(x), 1) .* (1 - tanh(Pe));
northBC = struct('A', northA, 'B', northB, 'C', northC);

eastA = @(y, t) y .* 0 + dirichletHackCoef;
eastB = @(y, t) y .* 0 + 1;
eastC = @(y, t) dirichletHackCoef .* ones(length(y), 1) .* (1 - tanh(Pe));
eastBC = struct('A', eastA, 'B', eastB, 'C', eastC);

southA = @(x, t) southA_problem3(x, t, dirichletHackCoef);
southB = @(x, t) x .* 0 + 1;
southC = @(x, t) southC_problem3(x, t, Pe, dirichletHackCoef);
southBC = struct('A', southA, 'B', southB, 'C', southC);

westA = @(y, t) y .* 0 + dirichletHackCoef;
westB = @(y, t) y .* 0 + 1;
westC = @(y, t) dirichletHackCoef .* ones(length(y), 1) .* (1 - tanh(Pe));
westBC = struct('A', westA, 'B', westB, 'C', westC);

% Construct initial condition
initialCondition = zeros(rows, columns);

% Initialise solver parameters
theta = 1;
advectionHandling = 'averaging';

storedTimeSteps = 1000;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'none', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = false;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
[tout_fine, yout_fine] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
