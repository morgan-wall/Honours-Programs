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
nodesX = 0:0.1:1;
nodesY = 0:0.2:2;
northBC = 0;
eastBC = 0;
southBC = 0;
westBC = 0;

initialCondition = zeros(length(nodesY), length(nodesX), 1);
initialCondition(5, 5) = 5;

storedTimeSteps = 1;

[tout, yout] = Solver(tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps);
