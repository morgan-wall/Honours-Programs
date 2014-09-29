function [root, iterations] = jacobian_free_newton_krylov(F, n, x0, ...
    max_iterations, rel_error_tol, gmres_error_tol, gmres_max_iter, ...
    restart_value, precond_type, omega)

%% Newton's Method

% initialise solver parameters
current_iteration = 0;
result = realmax;
x = x0;

while (current_iteration <= max_iterations && norm(result) > rel_error_tol)
    
    % solve the linear system using GMRES
    delta_x = jacobian_free_gmres(F, n, x, F(x), x0, gmres_max_iter, ...
        restart_value, gmres_error_tol, precond_type, omega);
    x = x - delta_x;
    
    result = F(x);
    current_iteration = current_iteration + 1;
end

root = x;
iterations = current_iteration;

end

%% Helper Functions

function [xk, H, g] = jacobian_free_gmres(F, n, u, b, x0, max_iter, ...
    restart_value, error_tol, type, omega)

% construct matrix preconditioner
diagonal_index = 0;
identity = eye(n);
rows = length(b);
Fu = F(u);
    
jacobian = sparse(zeros(n));
h = determine_newton_step_delta(u);
for j = 1:n
    delta_basis = identity(:, j);
    jacobian(:, j) = (F(u(:) + h .* delta_basis(:)) - Fu) ./ h;
end

switch(type)
    case 'ilu'
        setup.type = 'nofill';
        setup.milu = 'row';
        setup.droptol = 0.1;
        [L,U] = ilu(jacobian, setup);
    case 'jacobi'
        L = spdiags(diag(jacobian), diagonal_index, speye(rows));
        U = speye(length(b));
    case 'sor'
        L = spdiags(diag(jacobian) / omega, diagonal_index, tril(jacobian));
        U = speye(rows);
    case 'none'
        L = speye(rows);
        U = L;
end

% initialise solver paramaters
epsilon = sqrt(eps);

% initialise the current approximation of the solution
xk = x0;

% perform GMRES using the current approximation of the solution
i = 0;
residual = realmax;
while (residual > error_tol && i < max_iter)
    
    % initialise for the "first" Krylov subspace vector
    r0 = b - jacobian_vector_product(F, u, xk, epsilon);
    
    % generate the orthonormal basis
    [Q, H, g] = arnoldi_newton_krylov(F, u, L, U, r0, ...
        restart_value, error_tol, epsilon, n);
    
    % compute the minimiser (least squares problem)
    y = H(1:end-1, :) \ g(1:end-1);

    % update the solution
    x_precond = L * U * xk + Q * y;
    xk = U \ (L \ x_precond);
    
    % update loop parameters
    residual = g(end);
    i = i + 1;
end
end

function product = jacobian_vector_product(F, u, v, epsilon)
product = (F(u + epsilon * v) - F(u)) / epsilon;
end

function [delta] = determine_newton_step_delta(x)
if (norm(x) == 0)
    delta = sqrt(eps);
else
    delta = sqrt(eps) * norm(x);
end
end
