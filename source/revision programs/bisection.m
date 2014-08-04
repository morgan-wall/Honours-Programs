function [root, iterations] = bisection(f, neg_bound, pos_bound, ...
    max_iterations, abs_error_tol)

%% Error handling
root_value = 0;

if (f(neg_bound) >= root_value)
    throw(Mexception('bisection:invalidBound', ...
        'The negative bound provided is invalid.'));
end

if (f(pos_bound) <= root_value)
    throw(Mexception('bisection:invalidBound', ...
        'The positive bound provided is invalid.'));
end

%% Bisection method
current_iteration = 1;
result = realmax;

while (current_iteration <= max_iterations && abs(result) > abs_error_tol) 
    
    current_iteration = current_iteration + 1;
    
    mid_bound = (pos_bound - neg_bound) / 2 + neg_bound;
    result = f(mid_bound);
    root = mid_bound;
    
    % update the bounds
    if (result == root_value)
        break;
    elseif (result < root_value)
        neg_bound = mid_bound;
    elseif (result > root_value)
        pos_bound = mid_bound;
    end
end

iterations = current_iteration;

end