function [y] = sep_ROP_sparse(A,X,B)
    % Function to compute ROPs separated per batch; encoded as a sparse matrix
    % Args:
    %   A: (n,m,b) array
    %   X: (n^2*b*ncol) array
    %   B: (n,m,b) array
    % Returns:
    %   y: (m,b) array

    [n,m,b] = size(A);

    % reshape X to a matrix
    ncol = round(length(X)/(b*n*n));
    X = reshape(X, [b,n,n,ncol]);
    X = permute(X, [2,3,1,4]); % (b,n,n,ncol)
    X = reshape(X, [n^2*b, ncol]); % (n^2*b, ncol)

    % verify compatibility of dimensions
    if size(B,1) ~= n || size(B,2) ~= m || size(B,3) ~= b
        error('Dimension mismatch between A and B');
    end

    D = sep_ROP_sparse_mat(A,B); % (m*b,n^2*b)
    y = D*X; % (m*b,ncol)
    y = y(:); % (m*b*ncol)
end