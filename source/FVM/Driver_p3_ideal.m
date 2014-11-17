%% Driver_p3: produce the results for the SEB410 assignment.
% This script generates all the results included in Morgan Wall's report for
% the SEB410 assignment. The script contains solutions for the third problem
% aimed to test a vertex-centred finite volume method for solving generalised
% two-dimensional non-linear advection-diffusion equations.

clear all;
close all;

% Load the benchmark solution
load('p3_output_benchmark_fixed.mat');

benchmarkSolution_t7 = reshape(yout_fine(:, end), rows, columns);
size_benchmarkSolution_ty7 = size(benchmarkSolution_t7);
benchmarkSolution_t7(2:2:size_benchmarkSolution_ty7(1), :) = [];
benchmarkSolution_t7(:, 2:2:size_benchmarkSolution_ty7(2)) = [];
benchmarkSolution_coarseMesh = benchmarkSolution_t7;

% Initialise problem parameters
dt = 0.01;
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
xCount = 50;
xGeoParameters = struct('lowerIsGeometric', false, ...
    'upperIsGeometric', true, 'commonRatio', 0.9); 

yLower = 0;
yUpper = 1;
yCount = 25;
yGeoParameters = struct('lowerIsGeometric', false, ...
    'upperIsGeometric', true, 'commonRatio', 0.9);

[nodesX, nodesY] = GenerateStructuredMesh(xLower, xUpper, xCount, ...
    yLower, yUpper, yCount, xGeoParameters, yGeoParameters);
nodesY = flipud(nodesY);

rows = length(nodesY);
columns = length(nodesX);

% Construct finer mesh
dx = 0.01;
xNodesFine = xLower:dx:xUpper;
xNodesFine = sort(union(xNodesFine, nodesX), 'ascend');

collection = arrayfun(@(x) find(nodesX == x,1,'first'), xNodesFine , ...
    'UniformOutput', false);
xIndices = [];
for i = 1:length(collection)
    if (~isempty(collection{i}))
        xIndices = [xIndices i];
    end
end

dy = 0.01;
yNodesFine = yLower:dy:yUpper;
yNodesFine = sort(union(yNodesFine, nodesY), 'descend');

collection = arrayfun(@(x) find(nodesY == x,1,'first'), yNodesFine , ...
    'UniformOutput', false);
yIndices = [];
for i = 1:length(collection)
    if (~isempty(collection{i}))
        yIndices = [yIndices i];
    end
end

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
initialCondition_fine = zeros(length(yNodesFine), length(xNodesFine));

% Initialise solver parameters
theta = 1;
advectionHandling = 'averaging';

storedTimeSteps = 2 * ceil(tFinal / dt);

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-4, 'tolResidual', 1e-4);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 1000, ...
    'errorTol', 1e-4, 'preconditioningType', 'ilu', 'omega', 1.1);

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
[tout, yout, gmresIterations, nonlinearFnCalls, ...
    failed] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);

% Solve finer problem
dt = 0.001;

gmresParameters = struct('maxIterations', 1000, 'restartValue', 500, ...
    'errorTol', 1e-4, 'preconditioningType', 'ilu', 'omega', 1.1);

[tout_fine, yout_fine, gmresIterations_fine, nonlinearFnCalls_fine, ...
    failed_fine] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, xNodesFine, yNodesFine, northBC, eastBC, southBC, westBC, ...
    initialCondition_fine, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);

figure;

surf(nodesX, nodesY, reshape(yout(:, end), rows, columns), ...
    'EdgeColor','none','FaceColor', 'interp');
plotTitle = ['t = ' num2str(tout(end))];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');
colorbar;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls)]);
disp(['Failed: ' num2str(failed)]);

if (~failed)
    y = reshape(yout(:, end), rows, columns);
    y_fine = reshape(yout_fine(:, end), length(yNodesFine), length(xNodesFine));
    y_fine = y_fine(yIndices, xIndices);
    error = norm(y(:) - y_fine(:)) / sqrt(rows * columns);
    disp(['Error: ' num2str(error)]);
end