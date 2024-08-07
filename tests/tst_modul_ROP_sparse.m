% setup path
addpath('../lib/operators/ROP');

n = 3;
m = 5;
b = 2;
p = 4;
ncol = 3;

A_shape = struct();
A_shape.in = [n^2*b*ncol,1];
A_shape.out = [p*m*ncol,1];

% real
alpha = rand(n,m,b);
beta = rand(n,m,b);
Gamma = randsample([-1,1], p*b, true);
Gamma = reshape(Gamma, [p,b]);
A3 = @(X) modul_ROP3(alpha, X, beta, Gamma);
A_sparse = @(X) modul_ROP_sparse(alpha, X, beta, Gamma);

x = rand(A_shape.in);

assert(norm(A3(x)-A_sparse(x))/norm(A3(x))<1e-6, 'modul_ROP3 and modul_ROP_sparse are not equal');

% complex
alpha = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
beta = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
Gamma = randsample([-1,1], p*b, true);
Gamma = reshape(Gamma, [p,b]);
A3 = @(X) modul_ROP3(alpha, X, beta, Gamma);
A_sparse = @(X) modul_ROP_sparse(alpha, X, beta, Gamma);

x = rand(A_shape.in);

assert(norm(A3(x)-A_sparse(x))/norm(A3(x))<1e-6, 'modul_ROP3 and modul_ROP_sparse are not equal');