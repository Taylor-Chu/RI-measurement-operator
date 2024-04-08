% setup path
addpath('../lib/operators/ROP');

m = 5;
n = 3;

A_shape = struct();
A_shape.in = [n^2,1];
A_shape.out = [m,1];

alpha = rand(n,m);
beta = rand(n,m);

A = @(X) ROP(alpha, X, beta);
At = @(y) ROPt(alpha, y, beta);

adjoint_test(A, At, A_shape)