close all;
clear all;

%% Example 1
F = @(x) [x(1)^3 + x(2) - 1; x(2)^3 - x(1) + 1];
J = @(x) [3 * x(1)^2, 1 ; -1, 3 * x(2)^2];

x0 = [0.9; 0.1];

m = 3;

max_iterations = 100;
rel_error_tol = 10^-6;

[root, iterations] = ...
    shamanskii_system_exact(F, J, x0, m, max_iterations, rel_error_tol);

disp(['Iterations required: ', num2str(iterations), '.']);
