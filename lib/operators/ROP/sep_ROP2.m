function [y] = sep_ROP2(A,X,B)
    % Function to compute ROPs separated per batch
    % Args:
    %   A: (n,m,b) array
    %   X: (n^2*b) array
    %   B: (n,m,b) array
    % Returns:
    %   y: (m,b) array

    [n,m,b] = size(A);

    % reshape X to a matrix
    X = reshape(X, [b, n, n]);
    X = permute(X, [2,3,1]); % (n,n,b)

    % verify compatibility of dimensions
    if size(X,1) ~= n || size(X,2) ~= n || size(X,3) ~= b
        error('Dimension mismatch between A and X');
    end
    if size(B,1) ~= n || size(B,2) ~= m || size(B,3) ~= b
        error('Dimension mismatch between X and B');
    end

    XB = pagemtimes(X, B); % (n,n,b) x (n,m,b) = (n,m,b)

    y = sum(conj(A).*XB, 1); % (m,b)
    y = y/sqrt(m); % (m,b) normalization
    y = y(:); % (m*b)

end