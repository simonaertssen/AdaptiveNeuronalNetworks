clear all; close all; clc;
% In this script we will be simulating a simple fully connected network of
% theta neurons

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

%% Theta model parsameters:
F = @thetaneurons;
tnow = 0; tend = 10;
h = 0.01;

pars.N = 500;
pars.a_n = 0.6667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
seed = 1; rng(seed);
IC = randn(pars.N, 1)*0.1;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

%% Running the model:
[t, thetas] = DOPRI_threshold(F, tnow, tend, IC, h, pars);
[t, thetas] = spikesNaN(t, thetas);
plot(t, thetas)

