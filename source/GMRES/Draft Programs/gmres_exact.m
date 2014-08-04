function [x] = gmres_exact(A, b, x0)
%% Determine the solution to the linear system Ax = b using GMRES
% Determine the solution to the linear system Ax = b using the Generalised
% Minimum Residual Method.
%
%   Input:
%       A: a real nxn matrix.
%       b: a real nx1 vector.
%       x0: an approximation of the solution.
%
%   Output:
%       x: the solution vector.
%

%% Error handling
system_size = size(A);
rhs_size = length(b);
approx_size = length(x0);

rows = system_size(1);

if (rows ~= rhs_size || rhs_size ~= approx_size)
    throw(Mexception('gmres_exact:argumentDimensionsMismatch', ...
        'Input arguments have mismatched dimensions.'));
end

%% Generalised Minimum Residual Method (GMRES)

% determine the residual
r0 = b - A * x0;

% generate the orthonormal basis
[Q, H] = arnoldi(A, r0);

% solve the minimiser (least squares problem)
y = (H' * H) \ (H' * norm(r0) * eye(rows+1, 1));

% update the approximation to obtain the approximate solution
x = x0 + Q * y;

end