function [root, iterations] = newton(f, fdot, x0, ...
    max_iterations, abs_error_tol)

root_value = 0;
current_iteration = 1;
root = x0;
result = realmax;

while (current_iteration <= max_iterations && abs(result) > abs_error_tol)
    root = root - f(root) / fdot(root);
    
    result = f(root);
    if (result == root_value)
        break;
    end
    
    current_iteration = current_iteration + 1;
end

iterations = current_iteration;

end