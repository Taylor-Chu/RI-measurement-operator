% setup path
addpath('../lib/operators/ROP');

n = 20;
m = 40;
bs = logspace(0,5,20);

A_shape = struct();
A_shape.in = [n^2*b,1];
A_shape.out = [m*b,1];

times = zeros(length(bs),1);
times2 = zeros(length(bs),1);
times_adj = zeros(length(bs),1);

for i = 1:length(bs)
    i

    b = round(bs(i));

    alpha = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);
    beta = (rand(n,m,b)+1i*rand(n,m,b))/sqrt(2);

    A = @(X) sep_ROP(alpha, X, beta);
    A2 = @(X) sep_ROP2(alpha, X, beta);
    At = @(y) sep_ROPt(alpha, y, beta);

    X = randn(n^2*b,1)+1i*randn(n^2*b,1);

    tic;
    y = A(X);
    times(i) = toc;

    tic;
    y2 = A2(X);
    times2(i) = toc;

    tic;
    AtY = At(y);
    times_adj(i) = toc;
end

figure();
loglog(bs, times, 'r', bs, times2, 'b');
xlabel('b');
ylabel('time');
legend('sep ROP', 'sep ROP2');

figure();
loglog(bs, times_adj, 'r');
xlabel('b');
ylabel('time');
legend('sep ROP adjoint');