clear all;
close all;

tFinal = 0.1;
Dxx = @(phi) 0.1;
Dyy = @(phi) 0.1;
Vx = @(phi) 0;
Vy = @(phi) 0;
source = @(phi) 0;
theta = 0;
advectionHandling = 'averaging';
nodesX = 0:0.05:1;
nodesY = 1:-0.05:0;
northBC = 0;
eastBC = 0;
southBC = 0;
westBC = 0;

initialCondition = zeros(length(nodesY), length(nodesX), 1);
initialCondition(round(length(nodesY) / 2), round(length(nodesY) / 2)) = 5;

storedTimeSteps = 100;

[tout, yout] = Solver(tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps);

figure;
surf(yout(:, :, end));