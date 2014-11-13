function [B] = southB_problem3(x, t)
B = zeros(length(x), 1);
x_upper = x(x <= 0);
B(x <= 0) = 1;
B(x > 0) = 1;
end