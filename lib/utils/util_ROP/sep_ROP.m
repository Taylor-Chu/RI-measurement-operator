function [y] = sep_ROP(A,X,B)
    % Function to compute ROPs
    % Args:
    %   A: (n,m,b) array
    %   X: (n^2*b) array
    %   B: (n,m,b) array
    % Returns:
    %   y: (m,b) array

    [n,m,b] = size(A);

    % reshape X to a matrix
    X = reshape(X, [n^2, b]);

    % verify compatibility of dimensions
    if size(X,1) ~= n^2 || size(X,2) ~= b
        error('Dimension mismatch between A and X');
    end
    if size(B,1) ~= n || size(B,2) ~= m || size(B,3) ~= b
        error('Dimension mismatch between X and B');
    end

    y = zeros(m,b);
    for i = 1:b
        y(:,i) = ROP(A(:,:,i), X(:,i), B(:,:,i));
    end
    y = y(:);

end