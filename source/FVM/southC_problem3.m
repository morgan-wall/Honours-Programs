function [C] = southC_problem3(x, t, Pe, dirichletHackCoef)
C = zeros(length(x), 1);
x_upper = x(x <= 0);
C(x <= 0) = dirichletHackCoef .* (1 + tanh(Pe .* (2 .* x_upper + 1)));
C(x > 0) = 0;
end