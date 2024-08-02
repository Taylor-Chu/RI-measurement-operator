% setup path
addpath('../lib/operators/ROP');
addpath('../lib/utils');

m = 5;
n = 3;
b = 2;
p = 4;
ncol = 10;

A_shape = struct();
A_shape.in = [n^2*b*ncol,1];
A_shape.out = [p*m*ncol,1];

% real
alpha = rand(n,m,b);
beta = rand(n,m,b);
Gamma = randsample([-1,1], p*b, true);
Gamma = reshape(Gamma, [p,b]);
A = @(X) modul_ROP(alpha, X, beta, Gamma);
A3 = @(X) modul_ROP3(alpha, X, beta, Gamma);

x = rand(A_shape.in);
x = reshape(x, [n^2*b, ncol]);

Ax = zeros(m*p, ncol);
for i=1:ncol
    Ax(:,i) = A(x(:,i));
end
Ax = Ax(:);

x = x(:);

assert(norm(Ax-A3(x))/norm(A3(x))<1e-6, 'ROP2 and ROP3 are not equal');


% % complex
% alpha = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
% beta = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
% A = @(X) sep_ROP(alpha, X, beta);
% At = @(y) sep_ROPt(alpha, y, beta);
% adjoint_test(A, At, A_shape);