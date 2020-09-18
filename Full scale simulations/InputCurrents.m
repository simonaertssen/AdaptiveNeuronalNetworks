% In this script we will be studying the reponse of the theta neuron model
% on different types of input current.

clear all; close all; clc;

%% Excitability:
% The theta neuron model is of type 1 excitability, meaning that the frequency 
% of spikes goes up as the input current increases.
tnow = 0; tend = 10;
I = @(t) t;

pars.N = 1; pars.n = 2; pars.a_n = a_n(pars.n);
pars.eta0 = 0.5; pars.delta = 0.1; pars.K = -1;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

seed = 5; rng(seed); IC = -pi;
% pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
F = @thetaneuron; h = 0.001;

[t, x, tdr, xdr] = DOPRIthresh(F, 0, 10, IC, h, pars, pi);


[t, theta] = DOPRIthresh(@thetaneuron, 0, 10, IC, h, pars, pi);