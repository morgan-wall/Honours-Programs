close all;
clear all;

%% Example 1
F = @(x) [x(1)^3 + x(2) - 1; x(2)^3 - x(1) + 1];

x0 = [0.5; 0.5];
lambda_delta = 0.1;
abs_error_tol = 1e-12;

root = homotopy_simple(F, 2, x0, lambda_delta, abs_error_tol);

% output results
disp('Root: ');
disp(root);

disp('F(root): ');
disp(F(root));
