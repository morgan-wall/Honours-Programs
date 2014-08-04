function [root, iterations] = fixed_point(g, x0, max_iterations, ...
    abs_error_tol)

next_x = x0;
current_iteration = 0;
result = realmax;

while (current_iteration <= max_iterations && abs(result) > abs_error_tol)
    previous_x = next_x;
    next_x = g(previous_x);
    result = g(next_x) - next_x;
    
    current_iteration = current_iteration + 1;
end
    
root = next_x;
iterations = current_iteration;

end