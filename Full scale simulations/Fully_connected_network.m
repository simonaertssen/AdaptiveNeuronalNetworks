clear all; close all; clc;
% In this script we will be simulating a simple fully connected network of
% theta neurons

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

%% Theta model parameters:
F = @thetaneurons;
tnow = 0; tend = 10;
h = 0.0025;

pars.N = 100;
pars.a_n = 0.666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
seed = 1; rng(seed);
IC = randn(pars.N, 1)*2 + 1;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

%%
figure; hold on;
mu = 0; gamma = 0.5;
cvar = @(mu, gamma, N) tan(((rand(1,N) - mu)/gamma - 0.5)*pi);
h = cvar(mu, gamma, 1000);
histogram(h, 'Normalization','pdf');

xloc = mu + linspace(-2*gamma, 2*gamma, 100);
plot(xloc, cauchypdf(xloc, mu, gamma))

%% Running the model: theta and QIF models
[t, thetas] = DOPRI_threshold(F, tnow, tend, IC, h, pars);
drawthetas = spikesNaN(thetas);

z = orderparameter(thetas);

%% Drawing: the neurons over time
%ftheta = figure('Renderer', 'painters', 'Position', [100 800 1000 200]);

figure; box on;

yyaxis left
plot(t, drawthetas, '-', 'LineWidth', 1.5, 'Color', [0, 0, 1,  10/pars.N]);
ylim([-pi, pi])

ylabel('$\theta_i$','Interpreter','latex', 'FontSize', 20)

yyaxis right
plot(t, abs(z), '-k', 'LineWidth', 2.5);
ylim([0, 1])
ylabel('$\vert \bar{Z}(t) \vert$','Interpreter','latex', 'FontSize', 20)
xlabel('$t$','Interpreter','latex', 'FontSize', 20)
ax = gca; ax.YAxis(1).Color = [0, 0, 1, 1]; ax.YAxis(2).Color = 'k';


%% Drawing: the phase space
z = orderparameter(thetas);
f2 = figure('Renderer', 'painters', 'DefaultLegendFontSize',10);
plot(real(z), imag(z)) 




