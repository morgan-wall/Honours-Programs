function [x, iterations] = fixed_point_system(F, x0, ...
    max_iterations, abs_error_tol)

%% Error handling
system_size = length(F);
if (system_size ~= length(x0))
   throw(Mexception('fixed_point_system:invalidArgumentDimensions', ...
       'The arguments supplied have mismatched dimensions.'));
end

%% Fixed point iteration
previous_x = x0;
current_iteration = 0;
diff = realmax;

while (current_iteration <= max_iterations && abs(diff) > abs_error_tol)
   
    next_x = zeros(system_size, 1);
    for i = 1:system_size
       next_x(i) = F{i}(previous_x); 
    end
    
    diff = norm(next_x(:) - previous_x(:));
    previous_x = next_x;
    current_iteration = current_iteration + 1;
end

x = next_x;
iterations = current_iteration;

end