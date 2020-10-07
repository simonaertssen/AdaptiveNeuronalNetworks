clear all; close all; clc;

%% What random number generator is the fastest?
n = 100;
N = 10000;
netp = 0.2;
meandegree = netp*(N - 1);

%% Poisson:
tic;
Pp = @(x) poisspdf(x, meandegree)*N;
for i = 1:n
    degrees = poissrnd(meandegree, [N,1]);
end
elapsedtime = toc;
elapsedtime/n % = 0.0141s per sample

tic;
for i = 1:n
    ps = Pp(degrees);
end
elapsedtime = toc;
elapsedtime/n % = 0.0033s per sample

%% Binomial:
n = 2;
tic;
Pb = @(x) binopdf(x, N, netp)*N;
for i = 1:n
    degrees = binornd(N, netp, [N,1]);
end
elapsedtime = toc;
elapsedtime/n % = 1.3392s per sample

tic;
for i = 1:n
    ps = Pb(degrees);
end
elapsedtime = toc;
elapsedtime/n % = 0.0310s per sample

%% Numerical precision:
netp = 0.5;
for N = [100, 1000, 10000]
    degrees = poissrnd(meandegree, [N,1]);
    meandegree = netp*(N - 1);
    pps = Pp(degrees);
    pbs = Pb(degrees);
    immse(pps,pbs)
end
% Results: [45.8688, 0, 0]

%% Results:
% Poisson is the clear winner for speed and accurcay.