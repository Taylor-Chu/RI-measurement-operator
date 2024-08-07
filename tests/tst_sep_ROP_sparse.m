% setup path
addpath('../lib/operators/ROP');

n = 3;
m = 5;
b = 2;
ncol = 3;

A_shape = struct();
A_shape.in = [n^2*b*ncol,1];
A_shape.out = [m*b*ncol,1];

% real
alpha = rand(n,m,b);
beta = rand(n,m,b);
A3 = @(X) sep_ROP3(alpha, X, beta);
A_sparse = @(X) sep_ROP_sparse(alpha, X, beta);

x = rand(A_shape.in);

assert(norm(A3(x)-A_sparse(x))/norm(A3(x))<1e-6, 'sep_ROP3 and sep_ROP_sparse are not equal');

% complex
alpha = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
beta = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
A3 = @(X) sep_ROP3(alpha, X, beta);
A_sparse = @(X) sep_ROP_sparse(alpha, X, beta);

x = rand(A_shape.in);

assert(norm(A3(x)-A_sparse(x))/norm(A3(x))<1e-6, 'sep_ROP3 and sep_ROP_sparse are not equal');