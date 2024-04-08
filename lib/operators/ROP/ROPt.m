function [z] = ROPt(A,y,B)
    % Function to compute ROPs
    % Args:
    %   A: (n,m) array
    %   y: (m) array
    %   B: (n,m) array
    % Returns:
    %   z: (n^2) array

    [n,m] = size(A);

    % verify compatibility of dimensions
    if size(B,1) ~= n || size(B,2) ~= m
        error('Dimension mismatch between X and B');
    end
    if length(y) ~= m
        error('Dimension mismatch between y and A');
    end

    diagy = diag(y);

    z = A * diagy * B';
    z = z(:);
    z = z/sqrt(m);

end