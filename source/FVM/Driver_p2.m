%% Driver_p2: produce the results for the SEB410 assignment.
% This script generates all the results included in Morgan Wall's report for
% the SEB410 assignment. The script contains solutions for the second problem
% aimed to test a vertex-centred finite volume method for solving generalised
% two-dimensional non-linear advection-diffusion equations.

clear all;
close all;

% Initialise problem parameters
dt = 0.001;
tFinal = 10;

Dxx = @(phi, x, y, t) phi .^ 2;
Dyy = @(phi, x, y, t) 5 .* phi .^ 2;

C = 20;
Vx = @(phi, x, y, t) (C / 2) .* phi;
Vy = @(phi, x, y, t) (C / 2) .* phi;

source = @(phi, x, y, t) 100 * x .* y .* (1-x) .* (1-y) .* exp(2 * x.^4.5) ...
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
xCount = 50;
xGeoParameters = struct('lowerIsGeometric', false, ...
    'upperIsGeometric', false, 'commonRatio', 1); 

yLower = 0;
yUpper = 1;
yCount = 50;
yGeoParameters = struct('lowerIsGeometric', false, ...
    'upperIsGeometric', false, 'commonRatio', 1);

[nodesX, nodesY] = GenerateStructuredMesh(xLower, xUpper, xCount, ...
    yLower, yUpper, yCount, xGeoParameters, yGeoParameters);
nodesY = flipud(nodesY);

rows = length(nodesY);
columns = length(nodesX);

% Determine source term
[X, Y] = meshgrid(nodesX(:), nodesY(:));
sourceTerm = source(0, X(:), Y(:), 0);
% sourceTermMaple = sourceMaple(0, X(:), Y(:), 0);

% Initialise boundary conditions
dirichletHackCoef = 10000;

northA = @(x, t) x .* 0 + dirichletHackCoef;
northB = @(x, t) x .* 0 + 1;
northC = @(x, t) zeros(length(x), 1);
northBC = struct('A', northA, 'B', northB, 'C', northC);

eastA = @(y, t) y .* 0 + dirichletHackCoef;
eastB = @(y, t) y .* 0 + 1;
eastC = @(y, t) zeros(length(y), 1);
eastBC = struct('A', eastA, 'B', eastB, 'C', eastC);

southA = @(x, t) x .* 0 + dirichletHackCoef;
southB = @(x, t) x .* 0 + 1;
southC = @(x, t) zeros(length(x), 1);
southBC = struct('A', southA, 'B', southB, 'C', southC);

westA = @(y, t) y .* 0 + dirichletHackCoef;
westB = @(y, t) y .* 0 + 1;
westC = @(y, t) zeros(length(y), 1);
westBC = struct('A', westA, 'B', westB, 'C', westC);

% Construct initial condition
initialCondition = zeros(rows, columns);

% Construct steady state solution
steadyState = @(x, y) 10 * x .* y .* (1 - x) .* (1 - y) .* exp(x .^ 4.5);
steadyStateSolution = steadyState(X(:), Y(:));
analyticSolution = steadyStateSolution;

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

chordSteps = 3;

isGlobalised = false;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
[tout, yout] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);

% Output plots and metrics
figure;

for i = [1 2 3 10]
    tout(i)
    plot(nodesX, diag(flipud(reshape(yout(:, i), rows, columns))), 'LineWidth', 2);
    hold all;
end

plot(nodesX, diag(reshape(steadyStateSolution, rows, columns)), 'Xk', 'LineWidth', 2);

plotTitle = '';
title(plotTitle);
xlabel('x, y');
ylabel('Solution');

legend(['t = ' num2str(tout(1))], ['t = ' num2str(tout(2))], ...
    ['t = ' num2str(tout(3))], ['t = ' num2str(tout(10))], 'Steady State');

set(findall(gcf,'type','text'), 'fontSize', 12);
set(gca, 'fontSize', 11);

hold off;

figure;

surf(nodesX, nodesY, reshape(yout(:, 1), rows, columns));
plotTitle = 'Analytic Solution (Problem 2) where t = 0';
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

figure;

surf(nodesX, nodesY, reshape(sourceTerm, rows, columns));
plotTitle = 'Source Term (Problem 2)';
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');
colormap('parula');

set(findall(gcf,'type','text'), 'fontSize', 12);
set(gca, 'fontSize', 11);

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

%% Case: Backward-Euler with Averaging

disp('***** Begin: Backward-Euler (averaging) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'averaging';

storedTimeSteps = 250;

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
tic;
[tout_be_avg, yout_be_avg, gmresIterations_be_avg, nonlinearFnCalls_be_avg, ...
    failed_be_avg] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_avg)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_avg)]);
disp(['Failed: ' num2str(failed_be_avg)]);

if (~failed_be_avg)
    PlotSolutionAndAnalytic(nodesX, tout_be_avg, yout_be_avg, analyticSolution, ...
        rows, columns, 'be_avg');
end

disp('***** End: Backward-Euler (averaging) *****');

%% Case: Backward-Euler with upwinding

disp('***** Begin: Backward-Euler (upwinding) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

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
tic;
[tout_be_up, yout_be_up, gmresIterations_be_up, nonlinearFnCalls_be_up, ...
    failed_be_up] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_up)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_up)]);
disp(['Failed: ' num2str(failed_be_up)]);

if (~failed_be_up)
    PlotSolutionAndAnalytic(nodesX, tout_be_up, yout_be_up, analyticSolution, ...
        rows, columns, 'be_up');
end

disp('***** End: Backward-Euler (upwinding) *****');

%% Case: Crank-Nicolson with Averaging

disp('***** Begin: Crank-Nicolson (averaging) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'averaging';

storedTimeSteps = 250;

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
tic;
[tout_cn_avg, yout_cn_avg, gmresIterations_cn_avg, nonlinearFnCalls_cn_avg, ...
    failed_cn_avg] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_avg)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_avg)]);
disp(['Failed: ' num2str(failed_cn_avg)]);

if (~failed_cn_avg)
    PlotSolutionAndAnalytic(nodesX, tout_cn_avg, yout_cn_avg, analyticSolution, ...
        rows, columns, 'cn_avg');
end

disp('***** End: Crank-Nicolson (averaging) *****');

%% Case: Crank-Nicolson with upwinding

disp('***** Begin: Crank-Nicolson (upwinding) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

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
tic;
[tout_cn_up, yout_cn_up, gmresIterations_cn_up, nonlinearFnCalls_cn_up, ...
    failed_cn_up] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_up)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_up)]);
disp(['Failed: ' num2str(failed_cn_up)]);

if (~failed_cn_up)
    PlotSolutionAndAnalytic(nodesX, tout_cn_up, yout_cn_up, analyticSolution, ...
        rows, columns, 'cn_up');
end

disp('***** End: Crank-Nicolson (upwinding) *****');

%% Case: Backward-Euler with Averaging and backtracking

disp('***** Begin: Backward-Euler (averaging/backtracking) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'averaging';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'none', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_be_avg_linesearch, yout_be_avg_linesearch, ...
    gmresIterations_be_avg_linesearch, nonlinearFnCalls_be_avg_linesearch, ...
    failed_be_avg_linesearch] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_avg_linesearch)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_avg_linesearch)]);
disp(['Failed: ' num2str(failed_be_avg_linesearch)]);

if (~failed_be_avg_linesearch)
    PlotSolutionAndAnalytic(nodesX, tout_be_avg_linesearch, ...
        yout_be_avg_linesearch, analyticSolution, rows, columns, 'be_avg_linesearch');
end

disp('***** End: Backward-Euler (averaging/backtracking) *****');

%% Case: Backward-Euler with upwinding and backtracking

disp('***** Begin: Backward-Euler (upwinding/backtracking) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'none', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_be_up_linesearch, yout_be_up_linesearch, ...
    gmresIterations_be_up_linesearch, nonlinearFnCalls_be_up_linesearch, ...
    failed_be_up_linesearch] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc; 

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_up_linesearch)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_up_linesearch)]);
disp(['Failed: ' num2str(failed_be_up_linesearch)]);

if (~failed_be_up_linesearch)
    PlotSolutionAndAnalytic(nodesX, tout_be_up_linesearch, ...
        yout_be_up_linesearch, analyticSolution, rows, columns, 'be_up_linesearch');
end

disp('***** End: Backward-Euler (upwinding/backtracking) *****');

%% Case: Crank-Nicolson with averaging and backtracking

disp('***** Begin: Crank-Nicolson (averaging/backtracking) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'averaging';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'none', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_cn_avg_linesearch, yout_cn_avg_linesearch, ...
    gmresIterations_cn_avg_linesearch, nonlinearFnCalls_cn_avg_linesearch, ...
    failed_cn_avg_linesearch] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_avg_linesearch)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_avg_linesearch)]);
disp(['Failed: ' num2str(failed_cn_avg_linesearch)]);

if (~failed_cn_avg_linesearch)
    PlotSolutionAndAnalytic(nodesX, tout_cn_avg_linesearch, ...
        yout_cn_avg_linesearch, analyticSolution, rows, columns, 'cn_avg_linesearch');
end

disp('***** End: Crank-Nicolson (averaging/backtracking) *****');

%% Case: Crank-Nicolson with upwinding and backtracking

disp('***** Begin: Crank-Nicolson (upwinding/backtracking) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'none', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_cn_up_linesearch, yout_cn_up_linesearch, ...
    gmresIterations_cn_up_linesearch, nonlinearFnCalls_cn_up_linesearch, ...
    failed_cn_up_linesearch] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_up_linesearch)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_up_linesearch)]);
disp(['Failed: ' num2str(failed_cn_up_linesearch)]);

if (~failed_cn_up_linesearch)
    PlotSolutionAndAnalytic(nodesX, tout_cn_up_linesearch, ...
        yout_cn_up_linesearch, analyticSolution, rows, columns, 'cn_up_linesearch');
end

disp('***** End: Crank-Nicolson (upwinding/backtracking) *****');

%% Case: Backward-Euler with Averaging and Inexact Newton method (Assignment)

disp('***** Begin: Backward-Euler (averaging/assignment) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'averaging';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'assignment', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = false;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_be_avg_inexact_asgn, yout_be_avg_inexact_asgn, ...
    gmresIterations_be_avg_inexact_asgn, nonlinearFnCalls_be_avg_inexact_asgn, ...
    failed_be_avg_inexact_asgn] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_avg_inexact_asgn)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_avg_inexact_asgn)]);
disp(['Failed: ' num2str(failed_be_avg_inexact_asgn)]);

if (~failed_be_avg_inexact_asgn)
    PlotSolutionAndAnalytic(nodesX, tout_be_avg_inexact_asgn, ...
        yout_be_avg_inexact_asgn, analyticSolution, rows, columns, 'be_avg_inexact_asgn');
end

disp('***** End: Backward-Euler (averaging/assignment) *****');

%% Case: Backward-Euler with upwinding and Inexact Newton method (Assignment)

disp('***** Begin: Backward-Euler (upwinding/assignment) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'assignment', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = false;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_be_up_inexact_asgn, yout_be_up_inexact_asgn, ...
    gmresIterations_be_up_inexact_asgn, nonlinearFnCalls_be_up_inexact_asgn, ...
    failed_be_up_inexact_asgn] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_up_inexact_asgn)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_up_inexact_asgn)]);
disp(['Failed: ' num2str(failed_be_up_inexact_asgn)]);

if (~failed_be_up_inexact_asgn)
    PlotSolutionAndAnalytic(nodesX, tout_be_up_inexact_asgn, ...
        yout_be_up_inexact_asgn, analyticSolution, rows, columns, 'be_up_inexact_asgn');
end

disp('***** End: Backward-Euler (upwinding/assignment) *****');

%% Case: Crank-Nicolson with averaging and Inexact Newton method (Assignment)

disp('***** Begin: Crank-Nicolson (averaging/assignment) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'averaging';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'assignment', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = false;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_cn_avg_inexact_asgn, yout_cn_avg_inexact_asgn, ...
    gmresIterations_cn_avg_inexact_asgn, nonlinearFnCalls_cn_avg_inexact_asgn, ...
    failed_cn_avg_inexact_asgn] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_avg_inexact_asgn)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_avg_inexact_asgn)]);
disp(['Failed: ' num2str(failed_cn_avg_inexact_asgn)]);

if (~failed_cn_avg_inexact_asgn)
    PlotSolutionAndAnalytic(nodesX, tout_cn_avg_inexact_asgn, ...
        yout_cn_avg_inexact_asgn, analyticSolution, rows, columns, 'cn_avg_inexact_asgn');
end

disp('***** End: Crank-Nicolson (averaging/assignment) *****');

%% Case: Crank-Nicolson with upwinding and Inexact Newton method (Assignment)

disp('***** Begin: Crank-Nicolson (upwinding/assignment) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'assignment', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = false;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_cn_up_inexact_asgn, yout_cn_up_inexact_asgn, ...
    gmresIterations_cn_up_inexact_asgn, nonlinearFnCalls_cn_up_inexact_asgn, ...
    failed_cn_up_inexact_asgn] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_up_inexact_asgn)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_up_inexact_asgn)]);
disp(['Failed: ' num2str(failed_cn_up_inexact_asgn)]);

if (~failed_cn_up_inexact_asgn)
    PlotSolutionAndAnalytic(nodesX, tout_cn_up_inexact_asgn, ...
        yout_cn_up_inexact_asgn, analyticSolution, rows, columns, 'cn_up_inexact_asgn');
end

disp('***** End: Crank-Nicolson (upwinding/assignment) *****');

%% Case: Backward-Euler with Averaging and inexact Newton (choice 1)

disp('***** Begin: Backward-Euler (averaging/choice1) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'averaging';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice1', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = false;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_be_avg_inexact_c1, yout_be_avg_inexact_c1, ...
    gmresIterations_be_avg_inexact_c1, nonlinearFnCalls_be_avg_inexact_c1, ...
    failed_be_avg_inexact_c1] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_avg_inexact_c1)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_avg_inexact_c1)]);
disp(['Failed: ' num2str(failed_be_avg_inexact_c1)]);

if (~failed_be_avg_inexact_c1)
    PlotSolutionAndAnalytic(nodesX, tout_be_avg_inexact_c1, ...
        yout_be_avg_inexact_c1, analyticSolution, rows, columns, 'be_avg_inexact_c1');
end

disp('***** End: Backward-Euler (averaging/choice1) *****');

%% Case: Backward-Euler with upwinding and inexact Newton (choice 1)

disp('***** Begin: Backward-Euler (upwinding/choice1) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice1', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = false;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_be_up_inexact_c1, yout_be_up_inexact_c1, ...
    gmresIterations_be_up_inexact_c1, nonlinearFnCalls_be_up_inexact_c1, ...
    failed_be_up_inexact_c1] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_up_inexact_c1)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_up_inexact_c1)]);
disp(['Failed: ' num2str(failed_be_up_inexact_c1)]);

if (~failed_be_up_inexact_c1)
    PlotSolutionAndAnalytic(nodesX, tout_be_up_inexact_c1, ...
        yout_be_up_inexact_c1, analyticSolution, rows, columns, 'be_up_inexact_c1');
end

disp('***** End: Backward-Euler (upwinding/choice1) *****');

%% Case: Crank-Nicolson with averaging and inexact Newton (choice 1)

disp('***** Begin: Crank-Nicolson (averaging/choice1) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'averaging';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice1', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = false;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_cn_avg_inexact_c1, yout_cn_avg_inexact_c1, ...
    gmresIterations_cn_avg_inexact_c1, nonlinearFnCalls_cn_avg_inexact_c1, ...
    failed_cn_avg_inexact_c1] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_avg_inexact_c1)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_avg_inexact_c1)]);
disp(['Failed: ' num2str(failed_cn_avg_inexact_c1)]);

if (~failed_cn_avg_inexact_c1)
    PlotSolutionAndAnalytic(nodesX, tout_cn_avg_inexact_c1, ...
        yout_cn_avg_inexact_c1, analyticSolution, rows, columns, 'cn_avg_inexact_c1');
end

disp('***** End: Crank-Nicolson (averaging/choice1) *****');

%% Case: Crank-Nicolson with upwinding and inexact Newton (choice 1)

disp('***** Begin: Crank-Nicolson (upwinding/choice1) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice1', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = false;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_cn_up_inexact_c1, yout_cn_up_inexact_c1, ...
    gmresIterations_cn_up_inexact_c1, nonlinearFnCalls_cn_up_inexact_c1, ...
    failed_cn_up_inexact_c1] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_up_inexact_c1)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_up_inexact_c1)]);
disp(['Failed: ' num2str(failed_cn_up_inexact_c1)]);

if (~failed_cn_up_inexact_c1)
    PlotSolutionAndAnalytic(nodesX, tout_cn_up_inexact_c1, ...
        yout_cn_up_inexact_c1, analyticSolution, rows, columns, 'cn_up_inexact_c1');
end

disp('***** End: Crank-Nicolson (upwinding/choice1) *****');

%% Case: Backward-Euler with Averaging and inexact Newton (choice 2)

disp('***** Begin: Backward-Euler (averaging/choice2) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'averaging';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice2', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = false;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_be_avg_inexact_c2, yout_be_avg_inexact_c2, ...
    gmresIterations_be_avg_inexact_c2, nonlinearFnCalls_be_avg_inexact_c2, ...
    failed_be_avg_inexact_c2] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_avg_inexact_c2)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_avg_inexact_c2)]);
disp(['Failed: ' num2str(failed_be_avg_inexact_c2)]);

if (~failed_be_avg_inexact_c2)
    PlotSolutionAndAnalytic(nodesX, tout_be_avg_inexact_c2, ...
        yout_be_avg_inexact_c2, analyticSolution, rows, columns, 'be_avg_inexact_c2');
end

disp('***** End: Backward-Euler (averaging/choice2) *****');

%% Case: Backward-Euler with upwinding and inexact Newton (choice 2)

disp('***** Begin: Backward-Euler (upwinding/choice2) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice2', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = false;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_be_up_inexact_c2, yout_be_up_inexact_c2, ...
    gmresIterations_be_up_inexact_c2, nonlinearFnCalls_be_up_inexact_c2, ...
    failed_be_up_inexact_c2] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_up_inexact_c2)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_up_inexact_c2)]);
disp(['Failed: ' num2str(failed_be_up_inexact_c2)]);

if (~failed_be_up_inexact_c2)
    PlotSolutionAndAnalytic(nodesX, tout_be_up_inexact_c2, ...
        yout_be_up_inexact_c2, analyticSolution, rows, columns, 'be_up_inexact_c2');
end

disp('***** End: Backward-Euler (upwinding/choice2) *****');

%% Case: Crank-Nicolson with averaging and inexact Newton (choice 2)

disp('***** Begin: Crank-Nicolson (averaging/choice2) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'averaging';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice2', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = false;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_cn_avg_inexact_c2, yout_cn_avg_inexact_c2, ...
    gmresIterations_cn_avg_inexact_c2, nonlinearFnCalls_cn_avg_inexact_c2, ...
    failed_cn_avg_inexact_c2] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_avg_inexact_c2)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_avg_inexact_c2)]);
disp(['Failed: ' num2str(failed_cn_avg_inexact_c2)]);

if (~failed_cn_avg_inexact_c2)
    PlotSolutionAndAnalytic(nodesX, tout_cn_avg_inexact_c2, ...
        yout_cn_avg_inexact_c2, analyticSolution, rows, columns, 'cn_avg_inexact_c2');
end

disp('***** End: Crank-Nicolson (averaging/choice2) *****');

%% Case: Crank-Nicolson with upwinding and inexact Newton (choice 2)

disp('***** Begin: Crank-Nicolson (upwinding/choice2) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice2', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = false;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_cn_up_inexact_c2, yout_cn_up_inexact_c2, ...
    gmresIterations_cn_up_inexact_c2, nonlinearFnCalls_cn_up_inexact_c2, ...
    failed_cn_up_inexact_c2] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_up_inexact_c2)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_up_inexact_c2)]);
disp(['Failed: ' num2str(failed_cn_up_inexact_c2)]);

if (~failed_cn_up_inexact_c2)
    PlotSolutionAndAnalytic(nodesX, tout_cn_up_inexact_c2, ...
        yout_cn_up_inexact_c2, analyticSolution, rows, columns, 'cn_up_inexact_c2');
end

disp('***** End: Crank-Nicolson (upwinding/choice2) *****');

%% Case: Backward-Euler with Averaging and backtracking Inexact Newton method (Assignment)

disp('***** Begin: Backward-Euler (averaging/backtracking/assignment) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'averaging';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'assignment', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_be_avg_inexact_asgn_ls, yout_be_avg_inexact_asgn_ls, ...
    gmresIterations_be_avg_inexact_asgn_ls, nonlinearFnCalls_be_avg_inexact_asgn_ls, ...
    failed_be_avg_inexact_asgn_ls] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_avg_inexact_asgn_ls)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_avg_inexact_asgn_ls)]);
disp(['Failed: ' num2str(failed_be_avg_inexact_asgn_ls)]);

if (~failed_be_avg_inexact_asgn_ls)
    PlotSolutionAndAnalytic(nodesX, tout_be_avg_inexact_asgn_ls, ...
        yout_be_avg_inexact_asgn_ls, analyticSolution, rows, columns, 'be_avg_inexact_asgn_ls');
end

disp('***** End: Backward-Euler (averaging/backtracking/assignment) *****');

%% Case: Backward-Euler with upwinding and backtracking Inexact Newton method (Assignment)

disp('***** Begin: Backward-Euler (upwinding/backtracking/assignment) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'assignment', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_be_up_inexact_asgn_ls, yout_be_up_inexact_asgn_ls, ...
    gmresIterations_be_up_inexact_asgn_ls, nonlinearFnCalls_be_up_inexact_asgn_ls, ...
    failed_be_up_inexact_asgn_ls] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_up_inexact_asgn_ls)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_up_inexact_asgn_ls)]);
disp(['Failed: ' num2str(failed_be_up_inexact_asgn_ls)]);

if (~failed_be_up_inexact_asgn_ls)
    PlotSolutionAndAnalytic(nodesX, tout_be_up_inexact_asgn_ls, ...
        yout_be_up_inexact_asgn_ls, analyticSolution, rows, columns, 'be_up_inexact_asgn_ls');
end
    
disp('***** End: Backward-Euler (upwinding/backtracking/assignment) *****');

%% Case: Crank-Nicolson with averaging and backtracking Inexact Newton method (Assignment)

disp('***** Begin: Crank-Nicolson (averaging/backtracking/assignment) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'averaging';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'assignment', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_cn_avg_inexact_asgn_ls, yout_cn_avg_inexact_asgn_ls, ...
    gmresIterations_cn_avg_inexact_asgn_ls, nonlinearFnCalls_cn_avg_inexact_asgn_ls, ...
    failed_cn_avg_inexact_asgn_ls] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_avg_inexact_asgn_ls)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_avg_inexact_asgn_ls)]);
disp(['Failed: ' num2str(failed_cn_avg_inexact_asgn_ls)]);

if (~failed_cn_avg_inexact_asgn_ls)
    PlotSolutionAndAnalytic(nodesX, tout_cn_avg_inexact_asgn_ls, ...
        yout_cn_avg_inexact_asgn_ls, analyticSolution, rows, columns, 'cn_avg_inexact_asgn_ls');
end

disp('***** End: Crank-Nicolson (averaging/backtracking/assignment) *****');

%% Case: Crank-Nicolson with upwinding and backtracking Inexact Newton method (Assignment)

disp('***** Begin: Crank-Nicolson (upwinding/backtracking/assignment) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'assignment', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_cn_up_inexact_asgn_ls, yout_cn_up_inexact_asgn_ls, ...
    gmresIterations_cn_up_inexact_asgn_ls, nonlinearFnCalls_cn_up_inexact_asgn_ls, ...
    failed_cn_up_inexact_asgn_ls] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_up_inexact_asgn_ls)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_up_inexact_asgn_ls)]);
disp(['Failed: ' num2str(failed_cn_up_inexact_asgn_ls)]);

if (~failed_cn_up_inexact_asgn_ls)
    PlotSolutionAndAnalytic(nodesX, tout_cn_up_inexact_asgn_ls, ...
        yout_cn_up_inexact_asgn_ls, analyticSolution, rows, columns, 'cn_up_inexact_asgn_ls');
end
    
disp('***** End: Crank-Nicolson (upwinding/backtracking/assignment) *****');

%% Case: Backward-Euler with Averaging and backtracking inexact Newton (choice 1)

disp('***** Begin: Backward-Euler (averaging/backtracking/choice1) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'averaging';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice1', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_be_avg_inexact_c1_ls, yout_be_avg_inexact_c1_ls, ...
    gmresIterations_be_avg_inexact_c1_ls, nonlinearFnCalls_be_avg_inexact_c1_ls, ...
    failed_be_avg_inexact_c1_ls] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_avg_inexact_c1_ls)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_avg_inexact_c1_ls)]);
disp(['Failed: ' num2str(failed_be_avg_inexact_c1_ls)]);

if (~failed_be_avg_inexact_c1_ls)
    PlotSolutionAndAnalytic(nodesX, tout_be_avg_inexact_c1_ls, ...
        yout_be_avg_inexact_c1_ls, analyticSolution, rows, columns, 'be_avg_inexact_c1_ls');
end
    
disp('***** End: Backward-Euler (averaging/backtracking/choice1) *****');

%% Case: Backward-Euler with upwinding and backtracking inexact Newton (choice 1)

disp('***** Begin: Backward-Euler (upwinding/backtracking/choice1) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice1', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_be_up_inexact_c1_ls, yout_be_up_inexact_c1_ls, ...
    gmresIterations_be_up_inexact_c1_ls, nonlinearFnCalls_be_up_inexact_c1_ls, ...
    failed_be_up_inexact_c1_ls] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_up_inexact_c1_ls)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_up_inexact_c1_ls)]);
disp(['Failed: ' num2str(failed_be_up_inexact_c1_ls)]);

if (~failed_be_up_inexact_c1_ls)
    PlotSolutionAndAnalytic(nodesX, tout_be_up_inexact_c1_ls, ...
        yout_be_up_inexact_c1_ls, analyticSolution, rows, columns, 'be_up_inexact_c1_ls');
end
    
disp('***** End: Backward-Euler (upwinding/backtracking/choice1) *****');

%% Case: Crank-Nicolson with averaging and backtracking inexact Newton (choice 1)

disp('***** Begin: Crank-Nicolson (averaging/backtracking/choice1) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'averaging';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice1', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_cn_avg_inexact_c1_ls, yout_cn_avg_inexact_c1_ls, ...
    gmresIterations_cn_avg_inexact_c1_ls, nonlinearFnCalls_cn_avg_inexact_c1_ls, ...
    failed_cn_avg_inexact_c1_ls] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_avg_inexact_c1_ls)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_avg_inexact_c1_ls)]);
disp(['Failed: ' num2str(failed_cn_avg_inexact_c1_ls)]);

if (~failed_cn_avg_inexact_c1_ls)
    PlotSolutionAndAnalytic(nodesX, tout_cn_avg_inexact_c1_ls, ...
        yout_cn_avg_inexact_c1_ls, analyticSolution, rows, columns, 'cn_avg_inexact_c1_ls');
end
    
disp('***** End: Crank-Nicolson (averaging/backtracking/choice1) *****');

%% Case: Crank-Nicolson with upwinding and backtracking inexact Newton (choice 1)

disp('***** Begin: Crank-Nicolson (upwinding/backtracking/choice1) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice1', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_cn_up_inexact_c1_ls, yout_cn_up_inexact_c1_ls, ...
    gmresIterations_cn_up_inexact_c1_ls, nonlinearFnCalls_cn_up_inexact_c1_ls, ...
    failed_cn_up_inexact_c1_ls] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_up_inexact_c1_ls)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_up_inexact_c1_ls)]);
disp(['Failed: ' num2str(failed_cn_up_inexact_c1_ls)]);

if (~failed_cn_up_inexact_c1_ls)
    PlotSolutionAndAnalytic(nodesX, tout_cn_up_inexact_c1_ls, ...
        yout_cn_up_inexact_c1_ls, analyticSolution, rows, columns, 'cn_up_inexact_c1_ls');
end
    
disp('***** End: Crank-Nicolson (upwinding/backtracking/choice1) *****');

%% Case: Backward-Euler with Averaging and backtracking inexact Newton (choice 2)

disp('***** Begin: Backward-Euler (averaging/backtracking/choice2) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'averaging';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice2', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_be_avg_inexact_c2_ls, yout_be_avg_inexact_c2_ls, ...
    gmresIterations_be_avg_inexact_c2_ls, nonlinearFnCalls_be_avg_inexact_c2_ls, ...
    failed_be_avg_inexact_c2_ls] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_avg_inexact_c2_ls)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_avg_inexact_c2_ls)]);
disp(['Failed: ' num2str(failed_be_avg_inexact_c2_ls)]);

if (~failed_be_avg_inexact_c2_ls)
    PlotSolutionAndAnalytic(nodesX, tout_be_avg_inexact_c2_ls, ...
        yout_be_avg_inexact_c2_ls, analyticSolution, rows, columns, 'be_avg_inexact_c2_ls');
end
    
disp('***** End: Backward-Euler (averaging/backtracking/choice2) *****');

%% Case: Backward-Euler with upwinding and backtracking inexact Newton (choice 2)

disp('***** Begin: Backward-Euler (upwinding/backtracking/choice2) *****');

% Initialise solver parameters
theta = 1;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice2', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_be_up_inexact_c2_ls, yout_be_up_inexact_c2_ls, ...
    gmresIterations_be_up_inexact_c2_ls, nonlinearFnCalls_be_up_inexact_c2_ls, ...
    failed_be_up_inexact_c2_ls] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_be_up_inexact_c2_ls)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_be_up_inexact_c2_ls)]);
disp(['Failed: ' num2str(failed_be_up_inexact_c2_ls)]);

if (~failed_be_up_inexact_c2_ls)
    PlotSolutionAndAnalytic(nodesX, tout_be_up_inexact_c2_ls, ...
        yout_be_up_inexact_c2_ls, analyticSolution, rows, columns, 'be_up_inexact_c2_ls');
end
    
disp('***** End: Backward-Euler (upwinding/backtracking/choice2) *****');

%% Case: Crank-Nicolson with averaging and backtracking inexact Newton (choice 2)

disp('***** Begin: Crank-Nicolson (averaging/backtracking/choice2) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'averaging';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice2', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_cn_avg_inexact_c2_ls, yout_cn_avg_inexact_c2_ls, ...
    gmresIterations_cn_avg_inexact_c2_ls, nonlinearFnCalls_cn_avg_inexact_c2_ls, ...
    failed_cn_avg_inexact_c2_ls] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_avg_inexact_c2_ls)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_avg_inexact_c2_ls)]);
disp(['Failed: ' num2str(failed_cn_avg_inexact_c2_ls)]);

if (~failed_cn_avg_inexact_c2_ls)
    PlotSolutionAndAnalytic(nodesX, tout_cn_avg_inexact_c2_ls, ...
        yout_cn_avg_inexact_c2_ls, analyticSolution, rows, columns, 'cn_avg_inexact_c2_ls');
end
    
disp('***** End: Crank-Nicolson (averaging/backtracking/choice2) *****');

%% Case: Crank-Nicolson with upwinding and backtracking inexact Newton (choice 2)

disp('***** Begin: Crank-Nicolson (upwinding/backtracking/choice2) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'upwinding';

storedTimeSteps = 250;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'choice2', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

isGlobalised = true;
linesearchParam = 1e-4;
minLambda = 0.1;
maxLambda = 0.5;
maxBacktracks = 15;

% Solve problem
tic;
[tout_cn_up_inexact_c2_ls, yout_cn_up_inexact_c2_ls, ...
    gmresIterations_cn_up_inexact_c2_ls, nonlinearFnCalls_cn_up_inexact_c2_ls, ...
    failed_cn_up_inexact_c2_ls] = ...
    Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps, isGlobalised, ...
    linesearchParam, minLambda, maxLambda, maxBacktracks);
toc;

% Output plots and metrics
disp(['GMRES Iterations: ' num2str(gmresIterations_cn_up_inexact_c2_ls)]);
disp(['Nonlinear function calls: ' num2str(nonlinearFnCalls_cn_up_inexact_c2_ls)]);
disp(['Failed: ' num2str(failed_cn_up_inexact_c2_ls)]);

if (~failed_cn_up_inexact_c2_ls)
    PlotSolutionAndAnalytic(nodesX, tout_cn_up_inexact_c2_ls, ...
        yout_cn_up_inexact_c2_ls, analyticSolution, rows, columns, 'cn_up_inexact_c2_ls');
end

disp('***** End: Crank-Nicolson (upwinding/backtracking/choice2) *****');
