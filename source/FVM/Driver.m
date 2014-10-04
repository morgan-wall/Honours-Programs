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
advectionHandling = 'averaging';

newtonParameters = struct('maxIterations', 15, 'relErrorTol', 1e-10);

gmresParameters = struct('maxIterations', 200, 'restartValue', 100, ...
    'errorTol', 1e-10, 'preconditioningType', 'ilu', 'omega', 0);

%% Test: Gaussian Diffusion (G1) - Single Point Initial Condition

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
westBC = struct('A', 0, 'B', 1, 'C', 0);

% Construct initial condition
initialCondition = zeros(length(nodesY), length(nodesX));
initialCondition(round(length(nodesY) / 2), round(length(nodesY) / 2)) = 1;

% Solve problem
[tout, yout] = Solver(tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters);

% Output plots and metrics
figure;

surf(nodesX, nodesY, reshape(yout(:, end), rows, columns));
plotTitle = ['Test Problem (G1.1): Gaussian Diffusion (t = ' ...
    num2str(tout(end)) ')'];
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

totalMass = sum(sum(yout(2:end-1, end))) + (1/2) * yout(1, end) ...
    + (1/2) * yout(end, end);

disp(['Total at end time: ' num2str(totalMass) '.']);

figure;

surf(nodesX, nodesY, reshape(yout(:, 1), rows, columns));
plotTitle = 'Test Problem (G1.2): Gaussian Diffusion (t = 0)';
title(plotTitle);
xlabel('x');
ylabel('y');
zlabel('Solution');

%% Test: 1D Diffusion (H1) - Homogeneous Dirichlet Boundary Conditions

% Initialise temporal parameters
tFinal = 0.1;
storedTimeSteps = 100;

% Initialise equation parameters
Dxx = @(phi) 0;
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
northBC = struct('A', 1000, 'B', 1, 'C', 0);
eastBC = struct('A', 1000, 'B', 1, 'C', 0);
southBC = struct('A', 1000, 'B', 1, 'C', 0);
westBC = struct('A', 1000, 'B', 1, 'C', 0);

% Construct initial condition
initialCondition = zeros(length(nodesY), length(nodesX));
initialCondition(:, 1) = 1;

% Generate analytic solution
% Accountability note: this code was sourced from another student for 
% testing purposes only and should be rewritten or removed prior to
% submission.
D = 0.1;
xLength = 1;
analyticFn = @(n) (4 / (n * pi)) * sin(n * pi * nodesX ./ xLength) ...
    .* exp(-(n * pi / xLength)^2 * D * tFinal);
analyticSolution = zeros(1, length(nodesX));

for n = 1:2:31
    analyticSolution = analyticSolution + analyticFn(n);    
end

% Solve problem
[tout, yout] = Solver(tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters);

% Output plots and metrics
figure;

nodeCount = length(nodesY);
error = norm(yout(1:nodeCount, 2) - analyticSolution(:)) / sqrt(nodeCount);
disp(['H1.1 error: ' num2str(error)]);

plot(nodesX, analyticSolution, 'b');
hold on;
plot(nodesX, yout(1:nodeCount, 2), 'r*');
plotTitle = ['Test Problem (H1.1): 1D Diffusion (t = ' ...
    num2str(tout(end)) ')'];
title(plotTitle);
xlabel('x');
ylabel('y');
legend('Analytic Solution', 'Numeric Solution');

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
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters);

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
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters);

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
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters);

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
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters);

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
