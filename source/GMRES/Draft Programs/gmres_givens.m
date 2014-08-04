%% GMRESWithCGS: find solution to Ax = b using restarted GMRES.
% This approach uses k iterations of Arnoldi's method with C repetitions
% of GMRES.
function [xk] = gmres_givens(A, b, x0, k, C)

% initialise the current approximation of the solution
xk = x0;

% initialise orthonormal basis
m = max(size(A));
Q = zeros(m, k);
H = zeros(k+1, k);

% perform GMRES using the current approximation of the solution
for z = 1:C
    
    % initialise for the "first" Krylov subspace vector
    r0 = b - A * xk;
    beta = norm(r0);
    Q(:, 1) = r0 ./ beta;

    % generate the kth dimensional orthonormal basis for a subspace of the 
    % Krylov subspace
    for i = 1:k

        % obtain the "next" Krylov subspace vector
        w = A * Q(:, i);

        % project onto current orthonormal basis vectors
        for j = 1:i
            H(j, i) = Q(:, j)' * w;
            w = w - H(j, i) * Q(:, j);
        end

        % normalise and store the orthogonal basis vector
        % N.B. in obtaining H_k there is no need to store Q_{k+1}. Furthermore,
        % note that the norm of w is equal to zero if k was the dimension of A.
        % This is due to the fact that the previously generated orthonormal
        % basis vectors already describe the space in which x lies (i.e. R^k).
        % As such, there is no need to store the final row of H_k.
        if ((k == m && i < k) || (k ~= m))
            H(i+1, i) = norm(w);
            Q(:, i+1) = w / H(i+1, i);  
        end 
    end

    % compute the minimiser (least squares problem)
    y = (H' * H) \ (H' * beta * eye(k+1, 1));

    % update the solution
    xk = xk + Q(:, 1:k) * y;
end
end
