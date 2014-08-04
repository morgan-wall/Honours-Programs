function [root, iterations] = chord_system_exact(F, J, x0, ...
    max_iterations, rel_error_tol)

%% Error Handling
system_size = length(F);
[m, n] = size(J);
if (system_size ~= m)
   throw(Mexception('newton_system_exact:invalidArgumentDimensions', ...
       'The arguments have mismatched dimensions.'));
end

%% Newton's Method
current_iteration = 0;
result = realmax;
x = x0;

jacobian = J(x);

while (current_iteration <= max_iterations && norm(result) > rel_error_tol)
    delta_x = - jacobian \ F(x);
    x = x + delta_x;
    
    result = F(x);
    current_iteration = current_iteration + 1;
end

root = x;
iterations = current_iteration;

end