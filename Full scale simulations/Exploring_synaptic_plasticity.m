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
h = 0.005; tnow = h; tend = 100;

pars.N = 100;
pars.a_n = 0.666666666666666666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;

seed = 1; rng(seed);
IC = randn(pars.N,1);

pars.e = 0; %randcauchy(seed, pars.eta0, pars.delta, pars.N);

%% Results without synaptic scaling:
f_noSS = figure('Renderer', 'painters', 'Position', [50, 50, 800, 300]); hold on; box on;

winnames = ["Kempter1999Window", "Song2000Window", "ChrolCannon2012Window", "Waddington2014Window"];
colors = ["#0072BD", "#D95319", "#77AC30", "#A2142F"];
for i = 1:4
    name = winnames(i);
    plastopts = struct('SP', true, 'window', str2func(name), 'SS', false, 'IP', false);
    [t, thetas_full, K, Kmeans] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);
    drawthetas = spikesNaN(thetas_full);
    
    imrow(i) = subplot(2,2,i);
    xlim([tnow, tend]); 
    
    yyaxis left
    plot(t, drawthetas, '-', 'LineWidth', 1.5, 'Color', [0, 0, 1, 0.01], 'HandleVisibility', 'off')
    ylabel('$\theta_i$','Interpreter','latex', 'FontSize', labelfont)
    ylim([-pi, pi]); 

    yyaxis right
    plot(t, Kmeans, 'LineWidth', 2, 'Color', colors(i))
    ylabel('$\langle k \rangle$','Interpreter','latex', 'FontSize', labelfont)
    
    name = char(name);
    legend(sprintf('$$W(t)_%s$$', name(1)), 'Location', 'southwest', 'Interpreter', 'latex', 'FontSize', labelfont)

    ax = gca; ax.YAxis(1).Color = [0, 0, 1]; ax.YAxis(2).Color = 'k';
    xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
end

print(f_noSS, '../Figures/LearningWithoutScaling.png', '-dpng', '-r300')
close(f_noSS)






