clear all; close all; clc;
% In this script we will explore the synaptic plasticity of the theta
% neurons. In different plots we can observe different behaviour from the
% different learning windows, and we see we need synaptic scaling to
% account for an ever in- or decreasing mean synaptic strength.

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 15;
labelfont = 13;
export = true;

%% Make a GPU init handle:
if gpuDeviceCount > 0
    d = gpuDevice(gpuDeviceCount-1);
    disp(d)
end
initarray = make_GPUhandle();

%% Theta model parameters:
h = 0.005; tnow = h; tend = 80;

pars.N = 100;
pars.a_n = 0.666666666666666666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;

seed = 1; rng(seed);
IC = randn(pars.N,1);

pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

%% Test simulation:
clc
[t, thetas_full, K, Kmeans, info] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,@Waddington2014Window);
drawthetas = spikesNaN(thetas_full);

%% Plot results
figure; hold on; box on;
xlim([tnow, tend]); ylim([-pi, pi])
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)

yyaxis right
plot(t, Kmeans, '-k', 'LineWidth', 2)
ylabel('$\langle k \rangle$','Interpreter','latex', 'FontSize', labelfont)

yyaxis left
% plot(t, info, '-k', 'LineWidth', 0.1)
plot(t, drawthetas, '-', 'LineWidth', 1.5, 'Color', [0, 0, 1, 0.01])
ylabel('$\theta_i$','Interpreter','latex', 'FontSize', labelfont)

ax = gca; ax.YAxis(1).Color = [0, 0, 1]; ax.YAxis(2).Color = 'k';

