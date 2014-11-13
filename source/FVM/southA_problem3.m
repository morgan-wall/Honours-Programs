function [A] = southA_problem3(x, t, dirichletHackCoef)
A = zeros(length(x), 1);
A(x <= 0) = dirichletHackCoef;
A(x > 0) = 0;
end