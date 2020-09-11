clear all; close all;
tic;

%% Parameters
F = @thetaneurons;
tnow = 0; tend = 10;
h = 0.005;

pars.N = 10000;
IC = randn(pars.N, 1)*0.9 + 1;
pars.a_n = a_n(2);
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
seed = 0;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

%% Builtin Matlab functionality:
[t,x] = DOPRI45(F, tnow, tend, IC, h, pars);
plot(t,x)
toc;

% Elapsed time is 8.204944 seconds.
