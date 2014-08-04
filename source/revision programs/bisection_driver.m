close all;
clear all;

max_iterations = 100;
abs_error_tol = 10^-6;

%% Example 1
f = @(x) 0.9 * cos(x) - sqrt(x);

pos_bound = 0;
neg_bound = 1;

[root, iterations] = bisection(f, neg_bound, pos_bound, max_iterations, abs_error_tol);

disp(['A root exists at ', num2str(root), ' (found after ', num2str(iterations), ' iterations).']); 

%% Example 2
f = @(x) x^3 - 4 * x + 2;

neg_bound = 1;
pos_bound = 2;

[root, iterations] = bisection(f, neg_bound, pos_bound, max_iterations, abs_error_tol);

disp(['A root exists at ', num2str(root), ' (found after ', num2str(iterations), ' iterations).']); 
