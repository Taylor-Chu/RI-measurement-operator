% setup path
addpath('../lib/operators/ROP');

m = 4;
n = 3;
b = 5;

% Compare diagonal of dependent ROPs with separated ROPs
alpha = rand(n,m,b);
X = rand(n^2*b,1);
beta = rand(n,m,b);

A_sep = @(X) sep_ROP2(alpha, X, beta);
A_dep = @(X) dep_ROP(alpha, X, beta);

y_sep = A_sep(X);
y_dep = reshape(A_dep(X), [m,m,b]); % (m,m,1)
diag_y_dep = zeros(m,b); % (m,1)
for i=1:b
    diag_y_dep(:,i) = diag(y_dep(:,:,i));
end
diag_y_dep = diag_y_dep(:); % (m*b)

assert(norm(y_sep - diag_y_dep)/norm(y_sep) < 1e-10, 'Failed sep_ROP vs dep_ROP test');

% Adjoint test in real case
m = 5;
p = 6;
n = 3;
b = 2;

A_shape = struct();
A_shape.in = [n^2*b,1];
A_shape.out = [m*p*b,1];

% real
alpha = rand(n,m,b);
beta = rand(n,p,b);
A = @(X) dep_ROP(alpha, X, beta);
At = @(y) dep_ROPt(alpha, y, beta);
adjoint_test(A, At, A_shape);

% complex
alpha = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
beta = (rand(n,p,b)+1i*rand(n,p,b))/sqrt(2);
A = @(X) dep_ROP(alpha, X, beta);
At = @(y) dep_ROPt(alpha, y, beta);
adjoint_test(A, At, A_shape);