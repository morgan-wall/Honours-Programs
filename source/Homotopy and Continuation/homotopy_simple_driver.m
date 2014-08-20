close all;
clear all;

%% Example 1
F = @(x) [x(1)^3 + x(2) - 1; x(2)^3 - x(1) + 1];

x0 = [0.5; 0.5];
lambda_delta = 0.2;

root = homotopy_simple(F, 2, x0, lambda_delta);

% output results
disp('Root: ');
disp(root);

disp('F(root): ');
disp(F(root));
