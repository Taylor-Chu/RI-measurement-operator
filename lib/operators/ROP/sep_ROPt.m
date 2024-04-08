function [z] = sep_ROPt(A,y,B)
    % Function to compute ROPs
    % Args:
    %   A: (n,m,b) array
    %   y: (m*b) array
    %   B: (n,m,b) array
    % Returns:
    %   z: (n^2*b) array

    [n,m,b] = size(A);

    y = reshape(y,[m,b]);

    % verify compatibility of dimensions
    if size(B,1) ~= n || size(B,2) ~= m || size(B,3) ~= b
        error('Dimension mismatch between X and B');
    end
    if size (y,1) ~= m || size(y,2) ~= b
        error('Dimension mismatch between y and A');
    end

    z = zeros(b,n^2);
    for i=1:b 
        z(i,:) = ROPt(A(:,:,i),y(:,i),B(:,:,i));
    end
    z = z(:);

end