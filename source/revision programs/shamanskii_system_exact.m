function [root, iterations] = shamanskii_system_exact(F, J, x0, m, ...
    max_iterations, rel_error_tol)

%% Error Handling
system_size = length(F);
[n, k] = size(J);
if (system_size ~= n)
   throw(Mexception('newton_system_exact:invalidArgumentDimensions', ...
       'The arguments have mismatched dimensions.'));
end

%% Newton's Method
current_iteration = 0;
result = realmax;
x = x0;

while (current_iteration <= max_iterations && norm(result) > rel_error_tol)
    if (mod(current_iteration, m) == 0)
        jacobian = J(x);
    end
    
    delta_x = - jacobian \ F(x);
    x = x + delta_x;
    
    result = F(x);
    current_iteration = current_iteration + 1;
end

root = x;
iterations = current_iteration;

end