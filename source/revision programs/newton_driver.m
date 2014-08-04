clear all;
close all;

%% Example 1
f = @(x) x^3 - 4 * x + 2;
fdot = @(x) 3 * x^2 - 4;

x0 = 1.5;

max_iterations = 100;
abs_error_tol = 10^-6;

[root, iterations] = newton(f, fdot, x0, max_iterations, abs_error_tol);

disp(['A root exists at ', num2str(root), ' (found after ', num2str(iterations), ' iterations).']);
