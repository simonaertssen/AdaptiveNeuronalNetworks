%% Use a parfor loop to compute the function evaluation in parallel:
%% Setup:
N = 1000;
x = rand(N,1);
pars.a_n = matlab_a_n(2);
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
seed = 0;
pars.e = matlab_randcauchy(seed, pars.eta0, pars.delta, pars.N);

func = @(t, x) originalfunc(t, x, p.e, p.K/p.N, p.a_n);


K = zeros(
ButcherTableau = [
1/5         3/40    44/45   19372/6561      9017/3168       35/384
0           9/40    -56/15  -25360/2187     -355/33         0
0           0       32/9    64448/6561      46732/5247      500/1113
0           0       0       -212/729        49/176          125/192
0           0       0       0               -5103/18656     -2187/6784
0           0       0       0               0               11/84
0           0       0       0               0               0
];