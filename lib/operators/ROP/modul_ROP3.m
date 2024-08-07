function [z] = modul_ROP3(A,X,B,Gamma)
    % Function to compute modulated batched ROPs
    % Args:
    %   A: (n,m,b) array
    %   X: (n^2*b*ncol) array
    %   B: (n,m,b) array
    %   Gamma: (p,b) array
    %
    % Returns:
    %   z: (m,p) array

    [n,m,b] = size(A);
    ncol = round(length(X)/(n^2*b));

    % get the separated ROPs measurements
    y = sep_ROP3(A,X,B); % (m*b*ncol,)
    y = reshape(y, [m,b,ncol]); % (m,b,ncol)
    y = permute(y, [2,1,3]); % (b,m,ncol)
    y = reshape(y, [b,m*ncol]); % (b,m*ncol)

    % Modulate along the batches
    z = Gamma * y; % (p,m*ncol)
    z = z(:); % (p*m*ncol,)
    z = z/sqrt(b);

end