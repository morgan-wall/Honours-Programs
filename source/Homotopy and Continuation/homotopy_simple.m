function root = homotopy_simple(F, n, x0, lambda_delta, abs_error_tol)
%% homotopy_simple
% Attempt to find a root of a nonlinear system of equations using the
% embedding algorithm, which is a simple homotopy method.
%
%   Input:
%       F: the nonlinear system of equations (function handle).
%       n: the number of equations in the nonlinear system.
%       x0: an approximation of the solution (i.e. the zero point).
%       lambda_delta: the delta applied to the free variable in the
%           homotopy method.
%       abs_error_tol: the minimum absolute error imposed on the solution
%           of the Newton solver.
%
%   Output:
%       root: an approximation of the root for the nonlinear system.
%

% initialise problem parameters
min_lambda = 0;
max_lambda = 1;

% preallocate solution parameters
x = x0;

% initialise homotopy mapping
H = @(x, lambda) F(x) + (lambda - 1) * F(x0);

% initialise Newton solver parameters
max_iterations = 20;

%% Embedding Homotopy Method

for i = min_lambda:lambda_delta:max_lambda
    
    F_temp = @(x) H(x, i);
    
    [x, iterations] = newton_system_fde(F_temp, n, x, ... 
        max_iterations, abs_error_tol);
end

root = x;

end
