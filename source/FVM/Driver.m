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

% Construct analytic solution
xNodesAnalytic = linspace(xLower, xUpper, 50);
yNodesAnalytic = linspace(yLower, yUpper, 50);
rowsAnalytic = length(yNodesAnalytic);
columnsAnalytic = length(xNodesAnalytic);
[X, Y] = meshgrid(xNodesAnalytic, yNodesAnalytic);

analyticTimes = [0 0.25 0.5 0.75 1.0 1.25];
analyticSolutions = length(analyticTimes);
analyticSolution = zeros(rowsAnalytic * columnsAnalytic, analyticSolutions);

for i = 1:analyticSolutions
    analyticSolution(:, i) = phiAnalytic(X(:), Y(:), analyticTimes(i));
end

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

% Solve problem
[tout, yout] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps);

% Output plots and metrics
figure;

j = 1;
for i = 1:4
    
    subplot(2, 2, j);
    j = j + 1;
    
    surf(xNodesAnalytic, yNodesAnalytic, ...
        reshape(analyticSolution(:, i), rowsAnalytic, columnsAnalytic), ...
        'EdgeColor','none','FaceColor', 'interp');
    plotTitle = ['t = ' num2str(analyticTimes(i))];
    title(plotTitle);
    xlabel('x');
    ylabel('y');
    zlabel('Solution');
    
    zlim([0, 1]);
    
    set(findall(gcf,'type','text'), 'fontSize', 12);
    set(gca, 'fontSize', 11);
end

figure;

j = 1;
for i = 3:length(tout)
    
    subplot(2, 2, j);
    j = j + 1;
    
    plot(nodesX, diag(flipud( reshape(yout(:, i), rows, columns) )), 'LineWidth', 2);
    hold all;
    plot(xNodesAnalytic, diag(reshape(analyticSolution(:, i), ...
        rowsAnalytic, columnsAnalytic)), '--r', 'LineWidth', 2);
    
    plotTitle = ['t = ' num2str(tout(i))];
    title(plotTitle);
    xlabel('x, y');
    ylabel('Solution');

    ylim([0 0.4]);
    
    legend('Numeric', 'Analytic');

    set(findall(gcf,'type','text'), 'fontSize', 12);
        set(gca, 'fontSize', 11);
end

error = norm(yout(:, end) - analyticSolution(:, end)) / sqrt(length(yout(:, end)));

%% Problem 2: Non-linear convection-diffusion (with known steady-state)

% Initialise problem parameters
dt = 0.001;
tFinal = 15;

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

% Solve problem
[tout, yout] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps);

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

% figure;
% 
% surf(nodesX, nodesY, reshape(yout(:, 1), rows, columns));
% plotTitle = 'Analytic Solution (Problem 2) where t = 0';
% title(plotTitle);
% xlabel('x');
% ylabel('y');
% zlabel('Solution');
% 
% figure;
% 
% surf(nodesX, nodesY, reshape(sourceTerm, rows, columns));
% plotTitle = 'Source Term (Problem 2)';
% title(plotTitle);
% xlabel('x');
% ylabel('y');
% zlabel('Solution');
% colormap('parula');
% 
% set(findall(gcf,'type','text'), 'fontSize', 12);
% set(gca, 'fontSize', 11);
% 
% figure;
% 
% surf(nodesX, nodesY, reshape(yout(:, end), rows, columns));
% plotTitle = ['Numeric Solution (Problem 2) where t = ' num2str(tFinal)];
% title(plotTitle);
% xlabel('x');
% ylabel('y');
% zlabel('Solution');
% 
% figure;
% 
% surf(nodesX, nodesY, reshape(steadyStateSolution, rows, columns));
% plotTitle = 'Steady State Solution (Problem 2)';
% title(plotTitle);
% xlabel('x');
% ylabel('y');
% zlabel('Solution');

%% Problem 3: Non-linear convection-diffusion (incl. Peclet number)

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
advectionHandling = 'upwinding';

storedTimeSteps = 25;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'none', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

% Solve problem
[tout, yout] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps);

% Output plots and metrics
figure;
    
subplot(2, 2, 1);
surf(nodesX, nodesY, reshape(yout(:, 2), rows, columns), ...
    'EdgeColor','none','FaceColor', 'interp');
plotTitle = ['t = ' num2str(tout(2))];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

colormap('parula');

set(findall(gcf,'type','text'), 'fontSize', 12);
set(gca, 'fontSize', 11);

subplot(2, 2, 2);
surf(nodesX, nodesY, reshape(yout(:, 3), rows, columns), ...
    'EdgeColor','none','FaceColor', 'interp');
plotTitle = ['t = ' num2str(tout(3))];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

colormap('parula');

set(findall(gcf,'type','text'), 'fontSize', 12);
set(gca, 'fontSize', 11);

subplot(2, 2, 3);
surf(nodesX, nodesY, reshape(yout(:, 5), rows, columns), ...
    'EdgeColor','none','FaceColor', 'interp');
plotTitle = ['t = ' num2str(tout(5))];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

colormap('parula');

set(findall(gcf,'type','text'), 'fontSize', 12);
set(gca, 'fontSize', 11);

subplot(2, 2, 4);
surf(nodesX, nodesY, reshape(yout(:, end), rows, columns), ...
    'EdgeColor','none','FaceColor', 'interp');
plotTitle = ['t = ' num2str(tout(end))];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

colormap('parula');

set(findall(gcf,'type','text'), 'fontSize', 12);
set(gca, 'fontSize', 11);
