% Verify that the adjoint_test function is working correctly on a toy matrix example.
% Running this file should produce no output if the test passes.

n = 10;
m = 5;
A = randn(m,n);
A_shape = struct();
A_shape.in = [n,1];
A_shape.out = [m,1];
adjoint_test(A, A', A_shape);

% the following should throw an error
n = 10;
m = 10;
A = randn(m,n);
A_shape = struct();
A_shape.in = [n,1];
A_shape.out = [m,1];
adjoint_test(A, A, A_shape);