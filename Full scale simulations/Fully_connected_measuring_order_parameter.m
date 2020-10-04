clear all; close all; clc;
% In this script we will be testing the performance of different order
% parameters as suggested in Timme2017.

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 15;
labelfont = 18;

%% Theta model parameters:
F = @thetaneurons;
tnow = 0; tend = 5;
h = 0.01;

pars.N = 1000;
pars.a_n = 0.666667;
seed = 1; rng(seed);
IC = randn(pars.N, 1) + 1;

%% Network distributions:
A_fixeddegree = ones(pars.N, pars.N);



%% PSR state: one single stable node
pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

[t, thetas] = DOPRI_threshold(F, tnow, tend, IC, h, pars);

% Order parameters:
z = orderparameter(thetas);
degrees = sum(A_fixeddegree,2);
z_net = orderparameter_net(thetas, A_fixeddegree, sum(degrees));
z_mf = orderparameter_mf(thetas, degrees);

fspace = figure('Renderer', 'painters', 'Position', [50 800 800 200]);

hold on;
plot(t, abs(z), 'LineWidth', 2);
plot(t, abs(z_net), 'LineWidth', 2);
plot(t, abs(z_mf), 'LineWidth', 2);

xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
ylabel('$\| Z (t) \|$','Interpreter','latex', 'FontSize', labelfont)

legend('Kuramoto order parameter', 'Network order parameter', 'Mean field order parameter', 'FontSize', labelfont-5, 'Location', 'southeast')
removewhitspace();

