function [xk, H, g] = gmres_givens_short(A, b, x0, C)
%% Approximate the solution to the linear system Ax = b using GMRES
% Approximate the solution to the linear system Ax = b using C iterations
% of the Generalised Minimum Residual Method (GMRES) that generates k
% basis vectors using Arnoldi's method with Given's rotations.
%
%   Input:
%       A: a real nxn matrix.
%       b: a real nx1 vector.
%       x0: an approximation of the solution.
%       k: the number of vectors generated by the Arnoldi method each
%           iteration of the GMRES method.
%       C: the number of iterations of the GMRES method.
%
%   Output:
%       x: the solution vector.
%

% initialise the current approximation of the solution
xk = x0;

% perform GMRES using the current approximation of the solution
for z = 1:C
    
    % initialise for the "first" Krylov subspace vector
    r0 = b - A * xk;
    
    % generate the orthonormal basis
    [Q, H, g] = arnoldi_smart(A, r0);
    
    % compute the minimiser (least squares problem)
    y = H(1:end-1, :) \ g(1:end-1);

    % update the solution
    xk = xk + Q * y;
end
end