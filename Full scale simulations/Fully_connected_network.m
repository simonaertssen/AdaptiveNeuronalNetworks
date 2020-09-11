clear all; close all; clc;
% In this script we will be simulating a simple fully connected network of
% theta neurons

%% Setup
addpath('Functions');
addpath('Neuronmodel');
addpath('Adjacency matrices/');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')
rect = [700 400 500 250];

%% Parameters:
par.N = 10000;
par.n = 2;
par.a_n = a_n(par.n);
par.eta0 = 10.75; par.delta = 0.5; par.K = -9;
seed = 0; rng(seed);
IC = randn(par.N, 1)*2 + 1;
par.e = randcauchy(seed, par.eta0, par.delta, par.N);



%% Plots:
f1 = figure('Renderer', 'painters', 'Position', rect, 'DefaultLegendFontSize',10);
hold on 
%draw = plot(tdraw, xdraw(1,:), 'LineWidth', 1);
zplot = plot(t, abs(z), 'k--', 'LineWidth', 1.5);
Zplot = plot(T, abs(Z), 'b--','LineWidth', 1.5);
%zplot = plot(t, real(z), 'k--', 'LineWidth', 1.5);
%Zplot = plot(T, real(Z), 'b--','LineWidth', 1.5);
if isfield(par,'bif') == 1
    xline(par.bif.t, 'k-.','LineWidth', 1);
end
hold off
xlabel('Time $t$','Interpreter','latex')
ylabel('Phase $\theta$','Interpreter','latex')
legend([zplot, Zplot], 'order parameter $\| z \|$', 'order parameter $\| Z \|$','Interpreter','latex') % R2018b and later
close(f1)

f2 = figure('Renderer', 'painters', 'Position', [100, 500, 400, 400], 'DefaultLegendFontSize',10);
hold on 
zplot = plot(real(z), imag(z), 'k', 'LineWidth', 1.5);
Zplot = plot(real(Z), imag(Z), 'b', 'LineWidth', 1.5);
hold off
xlabel('Re{$z$}','Interpreter','latex')
ylabel('Im{$z$}','Interpreter','latex')

legend([zplot, Zplot], '$z$', '$Z$','Interpreter','latex') 

