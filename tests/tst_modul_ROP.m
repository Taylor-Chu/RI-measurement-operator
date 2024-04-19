% setup path
addpath('../lib/operators/ROP');
addpath('../lib/utils');

m = 5;
n = 3;
b = 2;
p = 4;

A_shape = struct();
A_shape.in = [n^2*b,1];
A_shape.out = [p*m,1];

% real
alpha = rand(n,m,b);
beta = rand(n,m,b);
Gamma = randsample([-1,1], p*b, true);
Gamma = reshape(Gamma, [p,b]);
A = @(X) modul_ROP(alpha, X, beta, Gamma);
At = @(y) modul_ROPt(alpha, y, beta, Gamma);
adjoint_test(A, At, A_shape);

% % complex
% alpha = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
% beta = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
% A = @(X) sep_ROP(alpha, X, beta);
% At = @(y) sep_ROPt(alpha, y, beta);
% adjoint_test(A, At, A_shape);