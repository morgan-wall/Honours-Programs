function [Q, H, g] = arnoldi_general(A, L, U, b, restart_value, error_tol)
%% Find orthonormal basis vectors for the Krylov subspace K(A, b)
% Determine orthonormal basis vectors for the Krylov subspace of b with 
% respect to the matrix A using Arnoldi's method for a specified restart
% value.
%
%   Input:
%       L: a real nxn matrix where elements above the diagonal are zero.
%       U: a real nxn matrix wher elements below the diagonal are zero.
%       b: a real nx1 vector.
%       restart_value: the restart value for the arnoldi process.
%       error_tol: the error tolerance imposed on the residual for
%           pre-emptively stopping the Arnoldi iteration within the GMRES
%           method.
%       type: the preconditioner used.
%
%   Output:
%       Q: a real nxk matrix containing the basis vectors.
%       H: the upper hessenberg matrix generated for determining the basis.
%       g: the rhs vector of the Given's relation which contains the
%           residual for the kth iteration of Arnoldi's algorithm.
%
%   Note: This function is used as part of GMRES
%

%% Error Handling
min_system_size = 2;
system_size = size(A);
rhs_size = length(b);

rows = system_size(1);
columns = system_size(2);

if (rows ~= rhs_size)
    throw(Mexception('arnoldi:argumentDimensionsMismatch', ...
        'A and b have mismatched dimensions.'));
elseif (columns < min_system_size)
    throw(Mexception('arnoldi:invalidArgument', ...
        ['A must have dimensions greater than ' num2str(min_system_size) '.']));
elseif (rows ~= columns)
    throw(Mexception('arnoldi:invalidArgument', 'A must be a square matrix.'));
end

%% Restarted Arnoldi's Method

% initialise solution variables
Q = zeros(rows, restart_value);
H = zeros(restart_value+1, restart_value);
g = ones(restart_value+1, 1); 
g(1) = norm(b);

% initialise Given's rotations variables
c = zeros(restart_value, 1);
s = zeros(restart_value, 1);

% generate the first basis vector
Q(:, 1) = b / norm(b);

% generate each basis vector in the Krylov subspace
for i = 1:restart_value
    
    % determine the next Krylov subspace basis vector
    q = A * (U \ (L \ Q(:, i)));
    
    % orthogonalise the basis vector
    for j = 1:i
        H(j, i) = Q(:, j)' * q;
        q = q - H(j, i) * Q(:, j);
    end
    
    % normalise the basis vector
    H(i+1, i) = norm(q);
    Q(:, i+1) = q / H(i+1, i);
    
    % apply Given's rotations
    for k = 1:i
        
        % construct rotations for new vector in Hessenberg matrix
        if (k == i)
            s(i) = H(i+1, i) / sqrt(H(i, i)^2 + H(i+1, i)^2);
            c(i) = H(i, i) / sqrt(H(i, i)^2 + H(i+1, i)^2);
        end
        
        H_diag_temp = c(k) * H(k, i) + s(k) * H(k+1, i);
        H(k+1, i) = -s(k) * H(k, i) + c(k) * H(k+1, i);
        H(k, i) = H_diag_temp;
    end
    
    % hack (avoid round-off error)
    H(i+1, i) = 0.0;
    
    g_temp = g(i);
    g(i) = c(i) * g_temp;
    g(i+1) = - s(i) * g_temp;
    
    if (abs(g(i+1)) < error_tol || i == restart_value)
        Q = Q(:, 1:i);
        H = H(1:i+1, 1:i);
        g = g(1:i+1);
        break;
    end
end
end