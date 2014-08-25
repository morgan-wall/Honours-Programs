function root = homotopy_simple(F, n, x0, lambda_delta)
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

%% Simple Homotopy Method

for i = min_lambda:lambda_delta:max_lambda
    
    identity = eye(n);
    jacobian = zeros(n);
    h = determine_jacobian_approx_delta(x);
    for j = 1:n
        delta_basis = identity(:, j);
        jacobian(:, j) = (H(x(:) + h .* delta_basis(:), i) - H(x, i)) ./ h;
    end
    
    delta_x = - jacobian \ H(x, i);
    x = x + delta_x;
end

root = x;

end

%% Helper Functions

function [delta] = determine_jacobian_approx_delta(x)
if (norm(x) == 0)
    delta = sqrt(eps);
else
    delta = sqrt(eps) * norm(x);
end
end