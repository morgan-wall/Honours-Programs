function root = homotopy_simple(F, n, x0, lambda_delta)

% initialise problem parameters
min_lambda = 0;
max_lambda = 1;

% preallocate solution parameters
x = x0;

%% Simple Homotopy Method

for i = min_lambda:lambda_delta:max_lambda
    
    identity = eye(n);
    jacobian = zeros(n);
    h = determine_jacobian_approx_delta(x);
    for j = 1:n
        delta_basis = identity(:, j);
        jacobian(:, j) = (F(x(:) + h .* delta_basis(:)) - F(x)) ./ h;
    end
    
    delta_x = - jacobian \ F(x);
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