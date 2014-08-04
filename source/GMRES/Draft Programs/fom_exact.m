function [x] = fom_exact(A, b, x0)
%% Determine the solution to the linear system Ax = b using FOM
% Determine the solution to the linear system Ax = b using Full 
% Orthogonalisation Method (FOM).
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
    throw(Mexception('fom:argumentDimensionsMismatch', ...
        'Input arguments have mismatched dimensions.'));
end

%% Full Orthogonalisation Method (FOM)

% determine the residual
r0 = b - A * x0;

% generate the orthonormal basis
[Q, H] = arnoldi(A, r0);

% update the approximation to obtain the approximate solution
y = H(1:end-1, :) \ (norm(r0) * eye(rows, 1));
x = x0 + Q * y;

end