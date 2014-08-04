clear all;
close all;

%% Example 1

g = @(x) (0.9 * cos(x))^2;

x0 = 0.5;

max_iterations = 100;
abs_error_tol = 10^-6;

[root, iterations] = fixed_point(g, x0, max_iterations, abs_error_tol);

disp(['A root exists at ', num2str(root), ' (found after ', num2str(iterations), ' iterations).']);

%% Example 2
% TODO: note that g(x) does not satisfy the fixed point theorem (fails to
% converge). We can overcome this issue by rearranging the original
% equation (N.B. currently, g(x) is too steep).

g = @(x) (1/4) * x^3 + (1/2);

x0 = 1.7;

max_iterations = 100;
abs_error_tol = 10^-6;

[root, iterations] = fixed_point(g, x0, max_iterations, abs_error_tol);

disp(['A root exists at ', num2str(root), ' (found after ', num2str(iterations), ' iterations).']);
