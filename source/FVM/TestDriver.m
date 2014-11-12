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
theta = 1;
advectionHandling = 'upwinding';

dt = 0.001;

newtonParameters = struct('rebuildJacobianIterations', 5, ...
    'maxIterations', 10, 'tolUpdate', 1e-8, 'tolResidual', 1e-8);

gmresParameters = struct('maxIterations', 1000, 'restartValue', 80, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

forcingTermParameters = struct('maxForcingTerm', 0.9, 'type', 'none', ...
    'gamma', 0.9, 'alpha', 2);

safeguardParameters = struct('threshold', 0.1);

chordSteps = newtonParameters.maxIterations + 1;

%% Test: Gaussian Diffusion (G1) - Single Point Initial Condition

% Initialise temporal parameters
tFinal = 0.4;
storedTimeSteps = 100;

% Initialise equation parameters
Dxx = @(phi) ones(length(phi), 1) * 0.1; %Dxx = @(phi) phi .* 0;
Dyy = @(phi) ones(length(phi), 1) * 0.1; %Dyy = @(phi) phi .* 0;
Vx = @(phi) phi .* 0; Vx = @(phi) ones(length(phi), 1) * 0.1;
Vy = @(phi) phi .* 0; Vy = @(phi) ones(length(phi), 1) * 0.1;
source = @(x, y) x .* 0;

% Initialise mesh parameters
xLower = 0;
xUpper = 1;
xCount = 85;
xGeoParameters = struct('lowerIsGeometric', false, 'upperIsGeometric', false, ...
    'commonRatio', 1); 

yLower = 0;
yUpper = 1;
yCount = 85;
yGeoParameters = struct('lowerIsGeometric', false, 'upperIsGeometric', false, ...
    'commonRatio', 1);

[nodesX, nodesY] = GenerateStructuredMesh(xLower, xUpper, xCount, ...
    yLower, yUpper, yCount, xGeoParameters, yGeoParameters);
nodesY = flipud(nodesY);

rows = length(nodesY);
columns = length(nodesX);

% Initialise boundary conditions
northC = @(x, t) x .* 0;
northBC = struct('A', -0.1, 'B', 0.1, 'C', northC);

eastC = @(y, t) y .* 0;
eastBC = struct('A', -0.1, 'B', 0.1, 'C', eastC);

southC = @(x, t) x .* 0;
southBC = struct('A', 0.1, 'B', 0.1, 'C', southC);

westC = @(y, t) y .* 0;
westBC = struct('A', 0.1, 'B', 0.1, 'C', westC);

% Construct initial condition
initialCondition = zeros(length(nodesY), length(nodesX));
initialCondition(round(length(nodesY) / 2), round(length(nodesX) / 2)) = 1;

% Solve problem
[tout, yout] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters, chordSteps);

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

% %% Test: 1D Diffusion (H1) - Homogeneous Dirichlet Boundary Conditions
% 
% % Initialise temporal parameters
% tFinal = 0.1;
% storedTimeSteps = 100;
% 
% % Initialise equation parameters
% Dxx = @(phi) phi .* 0;
% Dyy = @(phi) ones(length(phi), 1) .* 0.1;
% Vx = @(phi) phi .* 0;
% Vy = @(phi) phi .* 0;
% source = @(x, y) x .* 0;
% 
% % Initialise mesh parameters
% xLower = 0;
% xUpper = 1;
% xCount = 21;
% xGeoParameters = struct('lowerIsGeometric', true, 'upperIsGeometric', true, ...
%     'commonRatio', 1.1); 
% 
% yLower = 0;
% yUpper = 1;
% yCount = 21;
% yGeoParameters = struct('lowerIsGeometric', true, 'upperIsGeometric', true, ...
%     'commonRatio', 1.1);
% 
% [nodesX, nodesY] = GenerateStructuredMesh(xLower, xUpper, xCount, ...
%     yLower, yUpper, yCount, xGeoParameters, yGeoParameters);
% nodesY = flipud(nodesY);
% 
% nodesX = nodesX';
% nodesY = nodesY';
% 
% % Initialise boundary conditions
% northC = @(x, t) x .* 0;
% northBC = struct('A', 1000, 'B', 1, 'C', northC);
% 
% eastC = @(y, t) y .* 0;
% eastBC = struct('A', 1000, 'B', 1, 'C', eastC);
% 
% southC = @(x, t) x .* 0;
% southBC = struct('A', 1000, 'B', 1, 'C', southC);
% 
% westC = @(y, t) y .* 0;
% westBC = struct('A', 1000, 'B', 1, 'C', westC);
% 
% % Construct initial condition
% initialCondition = zeros(length(nodesY), length(nodesX));
% initialCondition(:, 1) = 1;
% 
% % Generate analytic solution
% % Accountability note: this code was sourced from another student for 
% % testing purposes only and should be rewritten or removed prior to
% % submission.
% D = 0.1;
% xLength = 1;
% analyticFn = @(n) (4 / (n * pi)) * sin(n * pi * nodesX ./ xLength) ...
%     .* exp(-(n * pi / xLength)^2 * D * tFinal);
% analyticSolution = zeros(1, length(nodesX));
% 
% for n = 1:2:31
%     analyticSolution = analyticSolution + analyticFn(n);    
% end
% 
% % Solve problem
% [tout, yout2] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
%     advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
%     initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
%     forcingTermParameters, safeguardParameters, chordSteps);
% 
% % Output plots and metrics
% figure;
% 
% nodeCount = length(nodesY);
% error = norm(yout2(1:nodeCount, 2) - analyticSolution(:)) / sqrt(nodeCount);
% disp(['H1.1 error: ' num2str(error)]);
% 
% plot(nodesX, analyticSolution, 'b');
% hold on;
% plot(nodesX, yout2(1:nodeCount, 2), 'r*');
% plotTitle = ['Test Problem (H1.1): 1D Diffusion (t = ' ...
%     num2str(tout(end)) ')'];
% title(plotTitle);
% xlabel('x');
% ylabel('y');
% legend('Analytic Solution', 'Numeric Solution');
% 
% %% Test: 1D Diffusion (H2) - Homogeneous Dirichlet Boundary Conditions
% 
% % Initialise temporal parameters
% tFinal = 0.1;
% storedTimeSteps = 100;
% 
% % Initialise equation parameters
% Dxx = @(phi) ones(length(phi), 1) .* 0.1;
% Dyy = @(phi) phi .* 0;
% Vx = @(phi) phi .* 0;
% Vy = @(phi) phi .* 0;
% source = @(x, y) x .* 0;
% 
% % Initialise mesh parameters
% xLower = 0;
% xUpper = 1;
% xCount = 21;
% xGeoParameters = struct('lowerIsGeometric', true, 'upperIsGeometric', true, ...
%     'commonRatio', 1.1); 
% 
% yLower = 0;
% yUpper = 1;
% yCount = 21;
% yGeoParameters = struct('lowerIsGeometric', true, 'upperIsGeometric', true, ...
%     'commonRatio', 1.1);
% 
% [nodesX, nodesY] = GenerateStructuredMesh(xLower, xUpper, xCount, ...
%     yLower, yUpper, yCount, xGeoParameters, yGeoParameters);
% nodesY = flipud(nodesY);
% 
% nodesX = nodesX';
% nodesY = nodesY';
% 
% rows = length(nodesY);
% columns = length(nodesX);
% 
% % Initialise boundary conditions
% northC = @(x, t) x .* 0;
% northBC = struct('A', 1000, 'B', 1, 'C', northC);
% 
% eastC = @(y, t) y .* 0;
% eastBC = struct('A', 1000, 'B', 1, 'C', eastC);
% 
% southC = @(x, t) x .* 0;
% southBC = struct('A', 1000, 'B', 1, 'C', southC);
% 
% westC = @(y, t) y .* 0;
% westBC = struct('A', 1000, 'B', 1, 'C', westC);
% 
% % Construct initial condition
% initialCondition = zeros(length(nodesY), length(nodesX));
% initialCondition(1, :) = 1;
% 
% % Generate analytic solution
% % Accountability note: this code was sourced from another student for 
% % testing purposes only and should be rewritten or removed prior to
% % submission.
% D = 0.1;
% xLength = 1;
% analyticFn = @(n) (4 / (n * pi)) * sin(n * pi * nodesX ./ xLength) ...
%     .* exp(-(n * pi / xLength)^2 * D * tFinal);
% analyticSolution = zeros(1, length(nodesX));
% 
% for n = 1:2:31
%     analyticSolution = analyticSolution + analyticFn(n);    
% end
% 
% % Solve problem
% [tout, yout] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
%     advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
%     initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
%     forcingTermParameters, safeguardParameters, chordSteps);
% 
% % Output plots and metrics
% figure;
% 
% nodeCount = length(nodesX);
% indices = 1:rows*columns;
% cIndices = indices - 1;
% northboundaryIndices = indices(mod(cIndices, rows) == 0);
% 
% error = norm(yout(northboundaryIndices, 2) - analyticSolution(:)) / sqrt(nodeCount);
% disp(['H2.1 error: ' num2str(error)]);
% 
% plot(nodesX, analyticSolution, 'b');
% hold on;
% plot(nodesX, yout(northboundaryIndices, 2), 'r*');
% plotTitle = ['Test Problem (H2.1): 1D Diffusion (t = ' ...
%     num2str(tout(end)) ')'];
% title(plotTitle);
% xlabel('x');
% ylabel('y');
% legend('Analytic Solution', 'Numeric Solution');
% 
% %% Test: Dirichlet & Neumann Boundary Conditions (N1) - Input at West face. 
% 
% % Initialise temporal parameters
% tFinal = 0.1;
% storedTimeSteps = 100;
% 
% % Initialise equation parameters
% Dxx = @(phi) 0.1;
% Dyy = @(phi) 0.1;
% Vx = @(phi) 0;
% Vy = @(phi) 0;
% source = @(x, y) x .* 0;
% 
% % Initialise mesh parameters
% nodesX = 0:0.05:1;
% nodesY = 1:-0.05:0;
% 
% rows = length(nodesY);
% columns = length(nodesX);
% 
% % Initialise boundary conditions
% northC = @(x, t) x .* 0;
% northBC = struct('A', 0, 'B', 1, 'C', northC);
% 
% eastC = @(y, t) y .* 0;
% eastBC = struct('A', 0, 'B', 1, 'C', eastC);
% 
% southC = @(x, t) x .* 0;
% southBC = struct('A', 0, 'B', 1, 'C', southC);
% 
% westC = @(y, t) y .* 0 + 1000;
% westBC = struct('A', 1000, 'B', 1, 'C', westC);
% 
% % Construct initial condition
% initialCondition = zeros(length(nodesY), length(nodesX));
% initialCondition(:, 1) = 1;
% 
% % Solve problem
% [tout, yout] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
%     advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
%     initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
%     forcingTermParameters, safeguardParameters, chordSteps);
% 
% % Output plots and metrics
% figure;
% 
% surf(nodesX, nodesY, reshape(yout(:, end), rows, columns));
% plotTitle = ['Test Problem (N1.1): Dirichlet & Neumann Boundary ' ...
%     'Conditions (t = ' num2str(tout(end)) ')'];
% title(plotTitle);
% xlabel('x');
% ylabel('y');
% zlabel('Solution');
% 
% figure;
% 
% surf(nodesX, nodesY, reshape(yout(:, 1), rows, columns));
% plotTitle = ['Test Problem (N1.2): Dirichlet & Neumann Boundary ' ...
%     'Conditions (t = 0 )'];
% title(plotTitle);
% xlabel('x');
% ylabel('y');
% zlabel('Solution');
% 
% %% Test: Dirichlet & Neumann Boundary Conditions (N2) - Input at North face. 
% 
% % Initialise temporal parameters
% tFinal = 0.1;
% storedTimeSteps = 100;
% 
% % Initialise equation parameters
% Dxx = @(phi) 0.1;
% Dyy = @(phi) 0.1;
% Vx = @(phi) 0;
% Vy = @(phi) 0;
% source = @(x, y) x .* 0;
% 
% % Initialise mesh parameters
% nodesX = 0:0.05:1;
% nodesY = 1:-0.05:0;
% 
% rows = length(nodesY);
% columns = length(nodesX);
% 
% % Initialise boundary conditions
% northC = @(x, t) x .* 0 + 1000;
% northBC = struct('A', 1000, 'B', 1, 'C', northC);
% 
% eastC = @(y, t) y .* 0;
% eastBC = struct('A', 0, 'B', 1, 'C', eastC);
% 
% southC = @(x, t) x .* 0;
% southBC = struct('A', 0, 'B', 1, 'C', southC);
% 
% westC = @(y, t) y .* 0;
% westBC = struct('A', 0, 'B', 1, 'C', westC);
% 
% % Construct initial condition
% initialCondition = zeros(length(nodesY), length(nodesX));
% initialCondition(1, :) = 1;
% 
% % Solve problem
% [tout, yout] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
%     advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
%     initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
%     forcingTermParameters, safeguardParameters, chordSteps);
% 
% % Output plots and metrics
% figure;
% 
% surf(nodesX, nodesY, reshape(yout(:, end), rows, columns));
% plotTitle = ['Test Problem (N2.1): Dirichlet & Neumann Boundary ' ...
%     'Conditions (t = ' num2str(tout(end)) ')'];
% title(plotTitle);
% xlabel('x');
% ylabel('y');
% zlabel('Solution');
% 
% figure;
% 
% surf(nodesX, nodesY, reshape(yout(:, 1), rows, columns));
% plotTitle = ['Test Problem (N2.2): Dirichlet & Neumann Boundary ' ...
%     'Conditions (t = 0 )'];
% title(plotTitle);
% xlabel('x');
% ylabel('y');
% zlabel('Solution');
% 
% %% Test: Dirichlet & Neumann Boundary Conditions (N3) - Input at East face. 
% 
% % Initialise temporal parameters
% tFinal = 0.1;
% storedTimeSteps = 100;
% 
% % Initialise equation parameters
% Dxx = @(phi) 0.1;
% Dyy = @(phi) 0.1;
% Vx = @(phi) 0;
% Vy = @(phi) 0;
% source = @(x, y) x .* 0;
% 
% % Initialise mesh parameters
% nodesX = 0:0.05:1;
% nodesY = 1:-0.05:0;
% 
% rows = length(nodesY);
% columns = length(nodesX);
% 
% % Initialise boundary conditions
% northC = @(x, t) x .* 0;
% northBC = struct('A', 0, 'B', 1, 'C', northC);
% 
% eastC = @(y, t) y .* 0 + 1000;
% eastBC = struct('A', 1000, 'B', 1, 'C', eastC);
% 
% southC = @(x, t) x .* 0;
% southBC = struct('A', 0, 'B', 1, 'C', southC);
% 
% westC = @(y, t) y .* 0;
% westBC = struct('A', 0, 'B', 1, 'C', westC);
% 
% % Construct initial condition
% initialCondition = zeros(length(nodesY), length(nodesX));
% initialCondition(:, end) = 1;
% 
% % Solve problem
% [tout, yout] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
%     advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
%     initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
%     forcingTermParameters, safeguardParameters, chordSteps);
% 
% % Output plots and metrics
% figure;
% 
% surf(nodesX, nodesY, reshape(yout(:, end), rows, columns));
% plotTitle = ['Test Problem (N3.1): Dirichlet & Neumann Boundary ' ...
%     'Conditions (t = ' num2str(tout(end)) ')'];
% title(plotTitle);
% xlabel('x');
% ylabel('y');
% zlabel('Solution');
% 
% figure;
% 
% surf(nodesX, nodesY, reshape(yout(:, 1), rows, columns));
% plotTitle = ['Test Problem (N3.2): Dirichlet & Neumann Boundary ' ...
%     'Conditions (t = 0 )'];
% title(plotTitle);
% xlabel('x');
% ylabel('y');
% zlabel('Solution');
% 
% %% Test: Dirichlet & Neumann Boundary Conditions (N4) - Input at South face. 
% 
% % Initialise temporal parameters
% tFinal = 0.1;
% storedTimeSteps = 100;
% 
% % Initialise equation parameters
% Dxx = @(phi) 0.1;
% Dyy = @(phi) 0.1;
% Vx = @(phi) 0;
% Vy = @(phi) 0;
% source = @(x, y) x .* 0;
% 
% % Initialise mesh parameters
% nodesX = 0:0.05:1;
% nodesY = 1:-0.05:0;
% 
% rows = length(nodesY);
% columns = length(nodesX);
% 
% % Initialise boundary conditions
% northC = @(x, t) x .* 0;
% northBC = struct('A', 0, 'B', 1, 'C', northC);
% 
% eastC = @(y, t) y .* 0;
% eastBC = struct('A', 0, 'B', 1, 'C', eastC);
% 
% southC = @(x, t) x .* 0 + 1000;
% southBC = struct('A', 1000, 'B', 1, 'C', southC);
% 
% westC = @(y, t) y .* 0;
% westBC = struct('A', 0, 'B', 1, 'C', westC);
% 
% % Construct initial condition
% initialCondition = zeros(length(nodesY), length(nodesX));
% initialCondition(end, :) = 1;
% 
% % Solve problem
% [tout, yout] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
%     advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
%     initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
%     forcingTermParameters, safeguardParameters, chordSteps);
% 
% % Output plots and metrics
% figure;
% 
% surf(nodesX, nodesY, reshape(yout(:, end), rows, columns));
% plotTitle = ['Test Problem (N4.1): Dirichlet & Neumann Boundary '...
%     'Conditions (t = ' num2str(tout(end)) ')'];
% title(plotTitle);
% xlabel('x');
% ylabel('y');
% zlabel('Solution');
% 
% figure;
% 
% surf(nodesX, nodesY, reshape(yout(:, 1), rows, columns));
% plotTitle = ['Test Problem (N4.2): Dirichlet & Neumann Boundary ' ...
%     'Conditions (t = 0 )'];
% title(plotTitle);
% xlabel('x');
% ylabel('y');
% zlabel('Solution');
