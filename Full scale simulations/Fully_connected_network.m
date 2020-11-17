clear all; close all; clc;
% In this script we will be simulating a simple fully connected network of
% theta neurons

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

%% Theta model parameters:
F = @thetaneurons;
tnow = 0; tend = 5;
h = 0.001;

pars.N = 10;
pars.a_n = 0.666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
seed = 1; rng(seed);
IC = wrapToPi(randn(pars.N, 1));
pars.e = 0.1;%randcauchy(seed, pars.eta0, pars.delta, pars.N);

%% Running the model: theta and QIF models
[t, thetas] = DOPRI_threshold(F, tnow, tend, IC, h, pars);

drawthetas = spikesNaN(thetas);

z = orderparameter(thetas);

%% Drawing: the neurons over time
if true
    ftheta = figure; box on;

    yyaxis left
    plot(t, drawthetas, '-', 'LineWidth', 1.5, 'Color', [0, 0, 1,  10/pars.N]);
    xlim([tnow, tend]); ylim([-pi, pi])

    ylabel('$\theta_i$','Interpreter','latex', 'FontSize', 20)

    yyaxis right
    plot(t, abs(z), '-k', 'LineWidth', 2.5);
    ylim([0, 1])
    ylabel('$\vert \bar{Z}(t) \vert$','Interpreter','latex', 'FontSize', 20)
    xlabel('$t$','Interpreter','latex', 'FontSize', 20)
    ax = gca; ax.YAxis(1).Color = [0, 0, 1, 1]; ax.YAxis(2).Color = 'k';
end

%% Drawing: the phase space
fspace = figure;
hold on;

scatter(real(z(:,1)), imag(z(:,1)), 50, 'k', 'filled', 'o', 'LineWidth',2);
plot(real(z), imag(z), '-k', 'LineWidth', 1.5);
phasespaceplot();




