function [root, iterations] = newton_system_fde(F, n, x0, ...
    max_iterations, rel_error_tol)

%% Newton's Method
m = length(x0);
current_iteration = 0;
result = realmax;
x = x0;

while (current_iteration <= max_iterations && norm(result) > rel_error_tol)
    
    identity = eye(n);
    jacobian = zeros(n);
    h = determine_newton_step_delta(x);
    for j = 1:n
        delta_basis = identity(:, j);
        jacobian(:, j) = (F(x(:) + h .* delta_basis(:)) - F(x)) ./ h;
    end
    
    delta_x = - jacobian \ F(x);
    x = x + delta_x;
    
    result = F(x);
    current_iteration = current_iteration + 1;
end

root = x;
iterations = current_iteration;

end

%% Helper Functions

function [delta] = determine_newton_step_delta(x)
if (norm(x) == 0)
    delta = sqrt(eps);
else
    delta = sqrt(eps) * norm(x);
end
end