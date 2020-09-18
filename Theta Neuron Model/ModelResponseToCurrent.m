% In this script we will be studying the reponse of the theta neuron model
% on different types of input current.

clear all; close all; clc;

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

%% Excitability:
% The theta neuron model is of type 1 excitability, meaning that the frequency 
% of spikes goes up as the input current increases.
tnow = 0; tend = 10;
inputcurrent = @(t) 10*t;
F = @thetaneuron; h = 0.001;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, -pi, h, inputcurrent);
drawthetas = spikesNaN(thetas);

fexcite = figure('Position', [50 800 1000 200]);
box on; hold on;

yyaxis left
ylim([-pi - 0.1, pi + 0.1]);
plot(t, thetas, ':k', 'LineWidth', 1)
plot(t, drawthetas, '-', 'LineWidth', 2.5, 'color', '#0072BD');
ylabel('$\vert \bar{Z}(t) \vert$','Interpreter','latex', 'FontSize', 20)

yyaxis right
plot(t, inputcurrent(t), '-', 'LineWidth', 1.5, 'color', '#A2142F');
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';
ylabel('$\vert \bar{Z}(t) \vert$','Interpreter','latex', 'FontSize', 20)
