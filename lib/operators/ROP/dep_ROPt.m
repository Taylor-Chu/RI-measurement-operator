function [AyB_star] = dep_ROPt(A,y,B)
    % Adjoint of dep_ROP
    % Args:
    %   A: (n,m,b) array
    %   y: (m*p*b) array
    %   B: (n,p,b) array
    % Returns:
    %   z: (b*n^2) array

    [n,m,b] = size(A);
    [~,p,~] = size(B);

    % verify compatibility of dimensions
    if size(B,1) ~= n || size(B,3) ~= b
        error('Dimension mismatch between X and B');
    end
    if size (y,1) ~= m*p*b
        error('Dimension mismatch between y and A');
    end

    y = reshape(y, [m,p,b]); % (m,p,b)
    B = permute(B, [2,1,3]); % (n,p,b) -> (p,n,b)

    yB_star = pagemtimes(y, conj(B)); % (m,n,b)
    AyB_star = pagemtimes(A, yB_star); % (n,n,b)
    AyB_star = permute(AyB_star, [3,1,2]); % (b,n,n)
    AyB_star = AyB_star/sqrt(m); % (n,n,b) normalization
    AyB_star = AyB_star(:); % (b*n^2)

end