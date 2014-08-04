close all;
clear all;

%% Example 1
F = @(x) [x(1)^3 + x(2) - 1; x(2)^3 - x(1) + 1];

x0 = [0.5; 0.5];

max_iterations = 100;
rel_error_tol = 10^-6;

[root, iterations] = ...
    newton_system_fde(F, 2, x0, max_iterations, rel_error_tol);

disp(['Iterations required: ', num2str(iterations), '.']);
