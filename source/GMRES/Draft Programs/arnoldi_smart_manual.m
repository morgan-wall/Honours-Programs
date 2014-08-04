function [Q, H, g] = arnoldi_smart_manual(A, b, restart_count)
%% Find orthonormal basis vectors for the Krylov subspace K(A, b)
% Determine the orthonormal basis for the Krylov subspace of b with 
% respect to the matrix A using Arnoldi's method with modified 
% Gram-Schmidt.
%
%   Input:
%       A: a real nxn matrix.
%       b: a real nx1 vector.
%
%   Output:
%       Q: a real nxk matrix containing the basis vectors.
%       H: the upper hessenberg matrix generated for determining the basis.
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

%% Arnoldi's Method

% set error bound on the norm of a basis vector
abs_error_tol = 10^-10;

% initialise solution variables
Q = zeros(rows);
H = zeros(rows+1, rows);
g = ones(rows+1, 1); 
g(1) = norm(b);

% initialise Given's rotations variables
c = zeros(rows, 1);
s = zeros(rows, 1);

% generate the first basis vector
Q(:, 1) = b / norm(b);

% generate each basis vector in the Krylov subspace
for i = 1:restart_count
    
    % determine the next Krylov subspace basis vector
    q = A * Q(:, i);
    
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
    
    if (abs(g(i+1)) < abs_error_tol)
        Q = Q(:, 1:i);
        H = H(1:i+1, 1:i);
        g = g(1:i+1);
        break;
    end
end
end