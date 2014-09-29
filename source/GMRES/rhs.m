function xdot = rhs(t, x, F_x0)
F_dot = @(x) [ 3 * x(1)^2, 1; -1, 3 * x(2)^2 ];
xdot = F_dot(x) * F_x0;
end