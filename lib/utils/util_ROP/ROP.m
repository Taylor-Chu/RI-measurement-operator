function [y] = ROP(A,X,B)
    % Function to compute ROPs
    % Args:
    %   A: (n,m) array
    %   X: (n^2) array
    %   B: (n,m) array
    % Returns:
    %   y: (m) array

    [n,m] = size(A);

    % reshape X to a matrix
    X = reshape(X, [n, n]);

    % verify compatibility of dimensions
    if size(X,1) ~= n || size(X,2) ~= n
        error('Dimension mismatch between A and X');
    end
    if size(B,1) ~= n || size(B,2) ~= m
        error('Dimension mismatch between X and B');
    end

    Rbeta_b = X * B; % (n,m)
    y = sum(conj(A) .* Rbeta_b, 1); % (m)
    y = transpose(y);
    y = y/sqrt(m);

end