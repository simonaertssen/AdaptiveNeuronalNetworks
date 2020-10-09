clear all; close all; clc;
% In this script we will be testing the performance of different order
% parameters as suggested in Timme2017.

%% Setup
addpath('../Functions');
addpath('../Mean Field Reductions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 15;
labelfont = 18;

%% Make a GPU init handle:
if gpuDeviceCount > 0
    d = gpuDevice(gpuDeviceCount);
    disp(d)
end
initarray = makeGPUinitHandle();

%% Theta model parameters:
F = @thetaneurons_adjacency;
tnow = 0; tend = 5;
h = 0.01;

pars.N = 1000;
pars.a_n = 0.666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
seed = 1; rng(seed);
IC = randn(pars.N, 1) + 1;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

%% Perform a full scale simulation of the dirac network:
diracpars = makeDiracPars(pars, 100);


%%
fspace = figure('Renderer', 'painters', 'Position', [50 800 800 200]);

hold on;
plot(t, abs(z), 'LineWidth', 2);
plot(t, abs(z_net), 'LineWidth', 2);
plot(t, abs(z_mf), 'LineWidth', 2);
plot(t, abs(z_link), 'LineWidth', 2);

xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
ylabel('$\| Z (t) \|$','Interpreter','latex', 'FontSize', labelfont)

legend('Kuramoto order parameter', 'Network order parameter', 'Mean field order parameter', 'Link field order parameter', 'FontSize', labelfont-5, 'Location', 'southeast')
removewhitspace();


