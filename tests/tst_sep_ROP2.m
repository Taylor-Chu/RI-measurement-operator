% setup path
addpath('../lib/operators/ROP');

n = 3;
m = 5;
b = 2;

A_shape = struct();
A_shape.in = [n^2*b,1];
A_shape.out = [m*b,1];

% real
alpha = rand(n,m,b);
beta = rand(n,m,b);
A = @(X) sep_ROP(alpha, X, beta);
A2 = @(X) sep_ROP2(alpha, X, beta);
At = @(y) sep_ROPt(alpha, y, beta);

x = rand(A_shape.in);
assert(norm(A(x)-A2(x))/norm(A(x))<1e-6, 'ROP and ROP2 are not equal');
adjoint_test(A, At, A_shape);

% complex
alpha = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
beta = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
A = @(X) sep_ROP(alpha, X, beta);
A2 = @(X) sep_ROP2(alpha, X, beta);
At = @(y) sep_ROPt(alpha, y, beta);

x = rand(A_shape.in);
assert(norm(A(x)-A2(x))/norm(A(x))<1e-6, 'ROP and ROP2 are not equal');
adjoint_test(A, At, A_shape);