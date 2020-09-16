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

pars.N = 100;
pars.a_n = 0.666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
seed = 1; rng(seed);
IC = randn(pars.N, 1)*2 + 1;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

%% Running the model: theta and QIF models
[t, thetas] = DOPRI_threshold(F, tnow, tend, IC, h, pars);
drawthetas = spikesNaN(thetas);

z = orderparameter(thetas);

%% Drawing: the neurons over time
ftheta = figure('Renderer', 'painters', 'Position', [100 800 1000 200]);

yyaxis left
plot(t, drawthetas, 'Color', [0, 0, 1, 0.05]);
ylabel('$\theta_i$','Interpreter','latex', 'FontSize', 20)
ylim([-pi, pi])

yyaxis right
plot(t, abs(z), '-k', 'LineWidth', 2.5);
ylabel('$\vert \bar{Z}(t) \vert$','Interpreter','latex', 'FontSize', 20)
ylim([0, 1])


%% Drawing: the phase space
z = orderparameter(thetas);
f2 = figure('Renderer', 'painters', 'DefaultLegendFontSize',10);
plot(real(z), imag(z)) 




