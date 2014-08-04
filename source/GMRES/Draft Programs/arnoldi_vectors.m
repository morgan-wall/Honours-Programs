function [Q, H] = arnoldi_vectors(A, b, k)
%% Find k orthonormal basis vectors for the Krylov subspace K(A, b)
% Determine the first k orthonormal basis for the Krylov subspace of b with 
% respect to the matrix A using Arnoldi's method with modified 
% Gram-Schmidt.
%
%   Input:
%       A: a real nxn matrix.
%       b: a real nx1 vector.
%       k: the number of basis vectors to generate.
%
%   Output:
%       Q: a real nxk matrix containing the basis vectors.
%       H: the upper hessenberg matrix generated for determining the basis
%           vectors.
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

if (k > rows)
   k = rows;
   disp(['Set k = ', num2str(rows), ...
       ' to avoid linear dependencies in basis vectors.']);
end

%% Arnoldi's Method

% set error bound on the norm of a basis vector
abs_error_tol = 10^-10;

Q = zeros(rows, k);
H = zeros(k+1, k);

% generate the first basis vector
Q(:, 1) = b / norm(b);

% generate each basis vector in the Krylov subspace
for i = 1:k
    
    % determine the next Krylov subspace basis vector
    q = A * Q(:, i);
    
    % orthogonalise the basis vector
    for j = 1:i
        H(j, i) = Q(:, j)' * q;
        q = q - H(j, i) * Q(:, j);
    end
    
    % normalise the basis vector
    q_norm = norm(q);
    if (q_norm > abs_error_tol)
        H(i+1, i) = q_norm;
        if (i ~= k)
            Q(:, i+1) = q / H(i+1, i);
        end
    else
        break;
    end
end

end