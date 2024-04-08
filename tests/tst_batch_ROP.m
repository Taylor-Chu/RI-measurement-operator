% setup path
addpath('../lib/operators/ROP');

m = 5;
n = 3;
b = 2;

A_shape = struct();
A_shape.in = [n^2*b,1];
A_shape.out = [m,1];

% real
alpha = rand(n,m,b);
beta = rand(n,m,b);
A = @(X) batch_ROP(alpha, X, beta);
At = @(y) batch_ROPt(alpha, y, beta);
adjoint_test(A, At, A_shape);

% complex
alpha = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
beta = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
A = @(X) batch_ROP(alpha, X, beta);
At = @(y) batch_ROPt(alpha, y, beta);
adjoint_test(A, At, A_shape);