function [y] = batch_ROP2(A,X,B)
    % Apply batched ROPs, which is separated ROPs summed over the time instants.
    % Args:
    %   A: (n,m,b) array
    %   X: (n^2*b) array
    %   B: (n,m,b) array
    % Returns:
    %   y: (m) array

    [~,m,b] = size(A);

    y = sep_ROP2(A,X,B); % (m*b)

    y = reshape(y,[m,b]); % (m,b)
    y = sum(y,2); % (m,1)
    y = y/sqrt(b); % normalize

end