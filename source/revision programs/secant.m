function [root, iterations] = secant(f, x0, x1, ...
    max_iterations, abs_error_tol)

root_value = 0;
current_iteration = 1;
previous_root = x0;
root = x1;
result = realmax;

while (current_iteration <= max_iterations && abs(result) > abs_error_tol)
    new_root = root ...
        - f(root) * (root - previous_root) / (f(root) - f(previous_root));
    previous_root = root;
    root = new_root;
    
    result = f(root);
    if (result == root_value)
        break;
    end
    
    current_iteration = current_iteration + 1;
end

iterations = current_iteration;

end