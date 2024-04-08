function [] = adjoint_test(A, A_star, A_shape)
    % Check if operator A_star is well the adjoint of operator A.

    % Args:
    %     A (torch.Tensor/LinearOperator): operator.
    %     A_star (torch.Tensor/LinearOperator): adjoint operator.
    %     A_shape (list): shape of the input and output of A. ex: [[2,3],[3,2]].

    if length(A_shape.in) > 2 || length(A_shape.out) > 2
        error('invalid dimension')
    end

    if A_shape.in(2) == 1 
        v = randn( A_shape.in );
    else 
        v = randn( [A_shape.in(1), A_shape.in(2)] );
    end

    if A_shape.out(2) == 1
        u = randn( A_shape.out );
    else
        u = randn( [A_shape.out(1), A_shape.out(2)] );
    end

    % if A is matrix, create function that multiply A and A_star with a matrix
    if isa(A, 'double')
        A_op = @(x) A*x;
        At_op = @(y) A_star*y;
    else 
        A_op = A;
        At_op = A_star;
    end

    score = abs(dot(u, A_op(v)) - dot(At_op(u), v))/(norm(u).*norm(A_op(v)));

    assert(score < 1e-8, sprintf('A and At are not adjoint, |<u,Av>-<A*u,v>|=%.2e||u|| ||Av||', score))

end