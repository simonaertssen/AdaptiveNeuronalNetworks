clear all; close all;
tic;

%% Parameters
F = @matlab_thetaneurons;
tnow = 0; tend = 10;
h = 0.005;

pars.N = 10000;
IC = randn(pars.N, 1)*0.9 + 1;
pars.a_n = matlab_a_n(2);
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
seed = 0;
pars.e = matlab_randcauchy(seed, pars.eta0, pars.delta, pars.N);

%% ODE solver results:
[t,x] = matlab_DOPRI(F, tnow, tend, IC, h, pars);
figure;
plot(t,x)
toc;

% Elapsed time is 8.204944 seconds (MacbookPro 2015, 2,9 GHz Intel Core i5).
% Elapsed time is 14.607684 seconds (DTU cluster)
% Elapsed time is 8.5326 seconds (DTU sxm2sh - X, matlab -nodisplay -nojvm < main_matlab_testspeed.m)


%% 
options = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
[t,x2] = ode45(@(t,x) F(t, x, pars.e, pars.K/pars.N, pars.a_n), [tnow, tend], IC, options);
figure;
plot(t,x2)
toc;

% Elapsed time is 23.531773 seconds (MacbookPro 2015, 2,9 GHz Intel Core i5).
