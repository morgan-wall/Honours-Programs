function [root, iterations] = newton_gmres_simple_line_search(F, n, x0, ...
    max_iterations, rel_error_tol, gmres_error_tol, gmres_max_iter, ...
    restart_value, precond_type, omega, sigma)

%% Newton's Method

% initialise solver parameters
identity = eye(n);
current_iteration = 0;
result = realmax;
x = x0;

while (current_iteration <= max_iterations && norm(result) > rel_error_tol)
    
    Fx = F(x);
    
    % determine finite difference approx of Jacobian
    jacobian = zeros(n);
    h = determine_newton_step_delta(x);
    for j = 1:n
        delta_basis = identity(:, j);
        jacobian(:, j) = (F(x(:) + h .* delta_basis(:)) - Fx) ./ h;
    end
    
    % solve the linear system using GMRES
    delta_x = gmres_general(jacobian, Fx, x0, gmres_max_iter, ...
        restart_value, gmres_error_tol, precond_type, omega);
    
    % perform a simple line search (using Armijo rule)
    lambda = 1;
    x_temp = x - lambda * delta_x;
    while (norm(F(x_temp)) >= (1 - sigma * lambda) * norm(Fx))
        lambda = lambda / 2;
        x_temp = x - lambda * delta_x;
    end
    x = x_temp;
    
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