function [AXB] = dep_ROP(A,X,B)
    % Function to compute dependent ROPs separated per batch
    % Args:
    %   A: (n,m,b) array
    %   X: (b*n^2) array
    %   B: (n,p,b) array
    % Returns:
    %   y: (m,p,b) array

    [n,m,b] = size(A);

    % reshape X to a matrix
    X = reshape(X, [b, n, n]);
    X = permute(X, [2,3,1]); % (n,n,b)

    % verify compatibility of dimensions
    if size(X,1) ~= n || size(X,2) ~= n || size(X,3) ~= b
        error('Dimension mismatch between A and X');
    end
    if size(B,1) ~= n || size(B,3) ~= b
        error('Dimension mismatch between X and B');
    end

    XB = pagemtimes(X, B); % (n,n,b) x (n,p,b) = (n,p,b)

    A = permute(A, [2,1,3]); % (m,n,b)
    
    AXB = pagemtimes(conj(A), XB); % (m,n,b) x (n,p,b) = (m,p,b)
    AXB = AXB/sqrt(m); % (m,p,b) normalization
    AXB = AXB(:); % (m*p*b) vectorization

end