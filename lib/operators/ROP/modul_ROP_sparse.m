function [z] = modul_ROP_sparse(A,X,B,Gamma)
    % Function to compute modulated batched ROPs
    % Args:
    %   A: (n,m,b) array
    %   X: (n^2*b*ncol) array
    %   B: (n,m,b) array
    %   Gamma: (p,b) array
    %
    % Returns:
    %   z: (m,p) array

    [n,m,b] = size(A);
    ncol = round(length(X)/(n^2*b));
    p = size(Gamma,1);

    % Block-diagonal matrix computing the separated ROPs
    D = sep_ROP_sparse_mat(A,B); % (m*b,n^2*b)

    % Sparse matrix computing the Bernoulli modulations
    M = modul_ROP_sparse_mat(Gamma, m); % (m*p, b*m)

    % reshape X to a matrix
    ncol = round(length(X)/(b*n*n));
    X = reshape(X, [b,n,n,ncol]);
    X = permute(X, [2,3,1,4]); % (b,n,n,ncol)
    X = reshape(X, [n^2*b, ncol]); % (n^2*b, ncol)

    % Modulate along the batches
    z = M * D * X; % (m*p,ncol)
    z = reshape(z, [m,p,ncol]); % (m,p,ncol)
    z = permute(z, [2,1,3]); % (p,m,ncol)
    z = z(:); % (p*m*ncol,)

end