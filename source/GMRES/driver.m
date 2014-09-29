clear all;
close all;

x0 = [0.5 0.5]';

lambda_span = [0 1];

F = @(x) [x(1)^3 + x(2) - 1; x(2)^3 - x(1) + 1];

F_x0 = F(x0);

rhs_fn = @(t, x) rhs(t, x, F_x0);

[lambda_sol, x_sol] = ode15s(rhs_fn, lambda_span, x0);

plot(lambda_sol, x_sol);

