function [y] = sep_ROP3(A,X,B)
    % Function to compute ROPs separated per batch
    % Args:
    %   A: (n,m,b) array
    %   X: (n^2*b*ncol) array
    %   B: (n,m,b) array
    % Returns:
    %   y: (m,b) array

    [n,m,b] = size(A);

    % reshape X to a matrix
    ncol = round(length(X)/(b*n*n));
    X = reshape(X, [b, n, n, ncol]);
    X = permute(X, [1,3,2,4]); % (b,n,n,ncol)
    X = reshape(X, [b,n,n*ncol]); % (b,n,n*ncol)
    X = permute(X, [3,2,1]); % (n*ncol,n,b)

    % verify compatibility of dimensions
    if size(X,1) ~= n*ncol || size(X,2) ~= n || size(X,3) ~= b
        error('Dimension mismatch between A and X');
    end
    if size(B,1) ~= n || size(B,2) ~= m || size(B,3) ~= b
        error('Dimension mismatch between X and B');
    end

    XB = pagemtimes(X, B); % (n*ncol,n,b) x (n,m,b) = (n*ncol,m,b)
    XB = reshape(XB, [n, ncol, m, b]);
    XB = permute(XB, [2,1,3,4]); % (ncol,n,m,b)
    % tmp = reshape(XB, [ncol, n*m*b]);

    % Repeat A for ncol times
    A = A(:); % (n*m*b)
    A = repmat(A, [1, ncol]); % (n*m*b, ncol)
    A = reshape(A, [n, m, b, ncol]);  
    A = permute(A, [4,1,2,3]); % (ncol,n,m,b)

    y = sum(conj(A).*XB, 2); % (ncol,m,b)
    % Squeeze y along second dimension
    y = reshape(y, [ncol, m, b]);
    y = permute(y, [2,3,1]); % (m,b,ncol)
    y = y/sqrt(m); % (m,b,ncol) normalization
    y = y(:); % (m*b*ncol)

end