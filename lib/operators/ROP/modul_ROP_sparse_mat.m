function [M] = modul_ROP_sparse_mat(Gamma, m)
    % Compute the sparse matrix applying the Bernoulli modulations.
    % Args:
    %   Gamma: (p,b) array
    %   m: int
    %
    % Returns:
    %   M: (m*p,b*m) array

    [p,b] = size(Gamma);

    % Sparse matrix computing the Bernoulli modulations
    M = kron(Gamma, speye(m)); % (m*p, b*m)
    M = M/sqrt(b); % normalization
end