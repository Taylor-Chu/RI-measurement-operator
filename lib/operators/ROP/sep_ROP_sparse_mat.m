function [D] = sep_ROP_sparse_mat(A,B)
    % Compute the block-diagonal matrix applying the separated ROPs.
    % Args:
    %   A: (n,m,b) array
    %   B: (n,m,b) array
    % Returns:
    %   D: (m*b,n^2*b) array

    [n,m,b] = size(A);

    % verify compatibility of dimensions
    if size(B,1) ~= n || size(B,2) ~= m || size(B,3) ~= b
        error('Dimension mismatch between A and B');
    end

    % Initialize the block-diagonal matrix computing the separated ROPs
    D = sparse(m*b, n^2*b);  % Using sparse for efficiency

    for k = 1:b      
        % Create (m x n^2) matrix by repeating A(:,:,k) n times
        Ak = repmat(reshape(A(:,:,k), [n, 1, m]), [1, n, 1]); % (n,n,m)
        Bk = repmat(reshape(B(:,:,k), [n, 1, m]), [1, n, 1]); % (n,n,m)
        
        % Create the block
        block = conj(Ak) .* permute(Bk, [2,1,3]); % (n,n,m)
        block = reshape(block, [n^2, m]); % (n^2,m)
        block = block.'; % (m,n^2)
        
        % Place the block in the correct position in C
        row_start = (k-1)*m + 1;
        row_end = k*m;
        col_start = (k-1)*n^2 + 1;
        col_end = k*n^2;
        
        D(row_start:row_end, col_start:col_end) = block;
    end

    D = D/sqrt(m); % normalization
end