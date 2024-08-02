% setup path
addpath('../lib/operators/ROP');

n = 3;
m = 5;
b = 2;
ncol = 10;

A_shape = struct();
A_shape.in = [n^2*b*ncol,1];
A_shape.out = [m*b*ncol,1];

% real
alpha = rand(n,m,b);
beta = rand(n,m,b);
A2 = @(X) sep_ROP2(alpha, X, beta);
A3 = @(X) sep_ROP3(alpha, X, beta);

x = rand(A_shape.in);
x = reshape(x, [n^2*b, ncol]);

A2x = zeros(m*b, ncol);
for i=1:ncol
    A2x(:,i) = A2(x(:,i));
end
A2x = A2x(:);

x = x(:);

assert(norm(A2x-A3(x))/norm(A3(x))<1e-6, 'ROP2 and ROP3 are not equal');

% complex
alpha = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
beta = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
A2 = @(X) sep_ROP2(alpha, X, beta);
A3 = @(X) sep_ROP3(alpha, X, beta);

x = rand(A_shape.in);
x = reshape(x, [n^2*b, ncol]);

A2x = zeros(m*b, ncol);
for i=1:ncol
    A2x(:,i) = A2(x(:,i));
end
A2x = A2x(:);

x = x(:);

assert(norm(A2x-A3(x))/norm(A3(x))<1e-6, 'ROP2 and ROP3 are not equal');