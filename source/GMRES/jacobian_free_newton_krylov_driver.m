close all;
clear all;

%% Example 1
F = @(x) [x(1, :) .^ 3 + x(2, :) - 1; x(2, :) .^ 3 - x(1, :) + 1];

x0 = [0.5; 0.5];

max_iterations = 100;
rel_error_tol = 10^-6;
gmres_max_iter = 10;
restart_value = 2;
gmres_error_tol = 10^-10; 
precond_type = 'ilu'; 
omega = 0;
sigma = 1e-4;

[root, iterations] = jacobian_free_newton_krylov(F, 2, x0, max_iterations, ...
    rel_error_tol, gmres_error_tol, gmres_max_iter, ...
    restart_value, precond_type, omega);

disp(['Iterations required: ', num2str(iterations), '.']);
