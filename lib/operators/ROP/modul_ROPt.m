function [AyB_star] = modul_ROPt(A,z,B,Gamma)
    % Adjoint of modul_ROP
    % Args:
    %   A: (n,m,b) array
    %   z: (p*m) array
    %   B: (n,m,b) array
    %   Gamma: (p,b) array
    %
    % Returns:
    %   AyB_star: (n^2*b) array

    [~,m,b] = size(A);
    p = size(Gamma,1);

    z = reshape(z,[p,m]); % (p,m)
    y = transpose(z) * Gamma; % (m,b)

    AyB_star = sep_ROPt(A,y,B); % (n^2*b)
    AyB_star = AyB_star/sqrt(b);

end