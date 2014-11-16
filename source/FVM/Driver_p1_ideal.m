%% Driver_p1: produce the results for the SEB410 assignment.
% This script generates all the results included in Morgan Wall's report for
% the SEB410 assignment. The script contains solutions for the first problem
% aimed to test a vertex-centred finite volume method for solving generalised
% two-dimensional non-linear advection-diffusion equations.

clear all;
close all;

% Initialise problem parameters
dt = 0.01;
tFinal = 1.25;

DXX = 0.01;
Dxx = @(phi, x, y, t) x .* 0 + DXX;

DYY = 0.01;
Dyy = @(phi, x, y, t) x .* 0 + DYY;

VX = 0.8;
Vx = @(phi, x, y, t) x .* 0 + VX;

VY = 0.8;
Vy = @(phi, x, y, t) x .* 0 + VY;

source = @(phi, x, y, t) x .* 0;

% Construct mesh
xLower = 0;
xUpper = 2;
xCount = 80;
xGeoParameters = struct('lowerIsGeometric', false, ...
    'upperIsGeometric', false, 'commonRatio', 1); 

yLower = 0;
yUpper = 2;
yCount = 80;
yGeoParameters = struct('lowerIsGeometric', false, ...
    'upperIsGeometric', false, 'commonRatio', 1);

[nodesX, nodesY] = GenerateStructuredMesh(xLower, xUpper, xCount, ...
    yLower, yUpper, yCount, xGeoParameters, yGeoParameters);
nodesY = flipud(nodesY);

rows = length(nodesY);
columns = length(nodesX);

% Initialise analytic solution
xC = 0.5;
yC = 0.5;

phiAnalytic = @(x, y, t) exp( -(x - VX * t - xC).^2 ./ (DXX * (4 * t + 1)) ...
    - (y - VY * t - yC).^2 ./ (DYY * (4 * t + 1)) ) ./(4 * t + 1);

% Initialise boundary conditions
dirichletHackCoef = 10000;

northA = @(x, t) x .* 0 + dirichletHackCoef;
northB = @(x, t) x .* 0 + 1;
northC = @(x, t) dirichletHackCoef .* phiAnalytic(x, yUpper, t);
northBC = struct('A', northA, 'B', northB, 'C', northC);

eastA = @(y, t) y .* 0 + dirichletHackCoef;
eastB = @(y, t) y .* 0 + 1;
eastC = @(y, t) dirichletHackCoef .* phiAnalytic(xUpper, y, t);
eastBC = struct('A', eastA, 'B', eastB, 'C', eastC);

southA = @(x, t) southA_problem3(x, t, dirichletHackCoef);
southB = @(x, t) x .* 0 + 1;
southC = @(x, t) dirichletHackCoef .* phiAnalytic(x, yLower, t);
southBC = struct('A', southA, 'B', southB, 'C', southC);

westA = @(y, t) y .* 0 + dirichletHackCoef;
westB = @(y, t) y .* 0 + 1;
westC = @(y, t) dirichletHackCoef .* phiAnalytic(xLower, y, t);
westBC = struct('A', westA, 'B', westB, 'C', westC);

% Construct initial condition
[X, Y] = meshgrid(nodesX(:), nodesY(:));
initialCondition = phiAnalytic(X(:), Y(:), 0);

% Construct analytic to mirror numeric solutions
[X, Y] = meshgrid(nodesX, nodesY);
storedTimeSteps = 25;
times = 0:storedTimeSteps*dt:tFinal;
analyticSolution = zeros(rows * columns, length(times));
for i = 1:length(times)
    analyticSolution(:, i) = phiAnalytic(X(:), Y(:), times(i));
end

phiAnalyticDiagSol = diag(flipud(reshape(analyticSolution(:, end), rows, columns)));

%% Case: Crank-Nicolson with Averaging (ideal)

disp('***** Begin: Crank-Nicolson (averaging) *****');

% Initialise solver parameters
theta = 1/2;
advectionHandling = 'averaging';

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-4, 'tolResidual', 1e-4);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 200, ...
    'errorTol', 1e-4, 'preconditioningType', 'ilu', 'omega', 0);

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

DetermineSecondaryErrorMetricsP1(yout_cn_avg, analyticSolution);

if (~failed_cn_avg)
    numericDiag = diag(flipud(reshape(yout_cn_avg(:, end), rows, columns)));
    [error] = CalculateErrorP1(numericDiag, phiAnalyticDiagSol);
    disp(['Error: ' num2str(error)]);
    
    figure;

    j = 1;
    for i = 3:length(tout_cn_avg)

        subplot(2, 2, j);
        j = j + 1;

        plot(nodesX, diag(flipud(reshape(yout_cn_avg(:, i), rows, columns) )), 'LineWidth', 2);
        hold all;
        plot(nodesX, diag(flipud(reshape(analyticSolution(:, i), rows, columns))), '--r', 'LineWidth', 2);

        plotTitle = ['t = ' num2str(tout_cn_avg(i))];
        title(plotTitle);
        xlabel('x, y');
        ylabel('Solution');

        ylim([0 0.4]);

        legend('Numeric', 'Analytic');

        set(findall(gcf,'type','text'), 'fontSize', 12);
            set(gca, 'fontSize', 11);
    end
end

disp('***** End: Crank-Nicolson (averaging) *****');
