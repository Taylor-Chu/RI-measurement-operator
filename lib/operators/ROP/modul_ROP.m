function [z] = modul_ROP(A,X,B,Gamma)
    % Function to compute modulated batched ROPs
    % Args:
    %   A: (n,m,b) array
    %   X: (n^2*b) array
    %   B: (n,m,b) array
    %   Gamma: (p,b) array
    %
    % Returns:
    %   z: (m,p) array

    [~,m,b] = size(A);

    % get the separated ROPs measurements
    y = sep_ROP(A,X,B); % (m*b,)
    y = reshape(y, [m,b]); % (m,b)

    % Modulate along the batches
    z = Gamma * transpose(y); % (p,m)
    z = z(:); % (p*m,)
    z = z/sqrt(b);

end