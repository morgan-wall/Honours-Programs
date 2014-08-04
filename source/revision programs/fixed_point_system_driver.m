clear all;
close all;

%% Example 1

f1 = @(x) x(2) / sqrt(5);
f2 = @(x) (1/4) * (sin(x(1)) + cos(x(2)));
F = { f1, f2 };

x0 = [ pi/2, pi/2 ];

max_iterations = 100;
abs_error_tol = 10^-6;

[root, iterations] = fixed_point_system(F, x0, max_iterations, abs_error_tol);

disp(['Iterations required: ', num2str(iterations)]);
