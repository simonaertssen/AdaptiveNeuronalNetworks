%clear all; close all; clc;
% In this script we will be simulating a simple fully connected network of
% theta neurons

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

%% Theta model parsameters:
F = @thetaneurons;
tnow = 0; tend = 10;
h = 0.0025;

pars.N = 8000;
pars.a_n = 0.6667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
seed = 1; rng(seed);
IC = randn(pars.N, 1)*2 + 1;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

%% Running the model:
[t, thetas] = DOPRI_threshold(F, tnow, tend, IC, h, pars);
drawthetas = spikesNaN(thetas);

% f1 = figure('Renderer', 'painters', 'DefaultLegendFontSize',10);
% plot(t, drawthetas)

%%
z = orderparameter(thetas);
f2 = figure('Renderer', 'painters', 'DefaultLegendFontSize',10);
plot(real(z), imag(z)) 




