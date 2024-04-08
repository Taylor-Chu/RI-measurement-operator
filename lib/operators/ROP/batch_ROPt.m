function [z] = batch_ROPt(A,y,B)
    % Function to compute ROPs
    % Args:
    %   A: (n,m,b) array
    %   y: (m) array
    %   B: (n,m,b) array
    % Returns:
    %   z: (n^2*b) array

    [n,m,b] = size(A);

    % verify compatibility of dimensions
    if size(B,1) ~= n || size(B,2) ~= m || size(B,3) ~= b
        error('Dimension mismatch between X and B');
    end
    if length(y) ~= m
        error('Dimension mismatch between y and A');
    end

    z = zeros(n^2,b);
    for i=1:b 
        z(:,i) = ROPt(A(:,:,i),y,B(:,:,i));
    end
    z = z(:);
    z = z/sqrt(b);

end