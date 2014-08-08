function [root, iterations] = newton_gmres(F, n, x0, ...
    max_iterations, rel_error_tol, gmres_error_tol, gmres_max_iter, ...
    restart_value, precond_type, omega)

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
    x = x - delta_x;
    
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