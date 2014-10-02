analytic = @(x, y, t) ( 5 / (4 * pi * t * sqrt(0.01 * 0.01)) ) ...
    * exp( (-1 / (4 * t)) * (x^2 / 0.01 + y^2 / 0.01) );

for i = 1:length(nodesX)
    for j = 1:length(nodesY)
        analytic_sol(i, j) = analytic(nodesX(i), nodesY(j), 0.1);
    end
end

figure;

surf(nodesX, nodesY, analytic_sol);
