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
h = 0.005; tnow = h; tend = 200;

pars.N = 100;
pars.a_n = 0.666666666666666666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;

seed = 2; rng(seed);
IC = randn(pars.N,1);

pars.e = 0; %randcauchy(seed, pars.eta0, pars.delta, pars.N);


%% Exploration:
STDP = struct('window', @Song2017Window, 'Kupdate', @(K, W) K + K.*W); 
plastopts = struct('SP', STDP, 'KMAX', 10);

K0 = initarray(10*ones(pars.N)); % For the learning windows

[t, thetas_full, K, Kmeans] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,K0,plastopts);
drawthetas = spikesNaN(thetas_full);
z = orderparameter(thetas_full);

figure; hold on; box on;
yyaxis left
plot(t, drawthetas, '-', 'LineWidth', 0.1, 'Color', [0, 0, 1, 1/pars.N], 'HandleVisibility', 'off')
plot(t, abs(z), '-k', 'LineWidth', 2)
ylabel('$\theta_i$','Interpreter','latex', 'FontSize', labelfont)
ylim([-pi, pi]);

yyaxis right
plot(t, Kmeans, 'LineWidth', 2, 'Color', "#0072BD")
ylabel('$\langle k \rangle$','Interpreter','latex', 'FontSize', labelfont)

ax = gca; ax.YAxis(1).Color = [0, 0, 1]; ax.YAxis(2).Color = "#0072BD";
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)



 %% Results without synaptic scaling:
f_noSS = figure('Renderer', 'painters', 'Position', [50, 50, 800, 300]); hold on; box on;

K0 = initarray(5*ones(pars.N)); % For the learning windows

winnames = ["Kempter1999Window", "Song2017Window", "ChrolCannon2012Window", "Waddington2014Window"];
colors = ["#0072BD", "#D95319", "#77AC30", "#A2142F"];
for i = 1:4
    name = winnames(i);
    
    STDP = struct('window', str2func(name), 'Kupdate', @(K, W) K + W);
    if name == "Kempter1999Window"
        STDP.w_i = 1.0e-5; STDP.w_o = - 1.0475*1.0e-5;
    end
    if name == "Song2017Window"
        STDP.Kupdate = @(K, W) K + K.*W;
    end
    
    plastopts = struct('SP', STDP, 'KMAX', 10);
    [t, thetas_full, ~, Kmeans] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,K0,plastopts);
    drawthetas = spikesNaN(thetas_full);
    z = orderparameter(thetas_full);

    imrow(i) = subplot(2,2,i);
    xlim([tnow, tend]); 
    
    yyaxis left
    rasterplot(t, drawthetas, labelfont);

%     plot(t, drawthetas, '-', 'LineWidth', 1.5, 'Color', [0, 0, 1, 0.01], 'HandleVisibility', 'off')
%     plot(t, abs(z), '-k', 'LineWidth', 2)
%     ylabel('$\theta_i$','Interpreter','latex', 'FontSize', labelfont)
%     ylim([-pi, pi]); 
%     set(gca,'YTick',-pi:pi/2:pi) 
%     set(gca,'YTickLabel',{'$$-\pi$$','$$-\frac{\pi}{2}$$','$$0$$','$$\frac{\pi}{2}$$','$$\pi$$'}, 'TickLabelInterpreter', 'latex')

    yyaxis right
    plot(t, Kmeans(1,:), 'LineWidth', 2, 'Color', colors(i))
    ylabel('$\langle k \rangle$','Interpreter','latex', 'FontSize', labelfont)
    
    name = char(name);
    title(sprintf('$$W(t)_%s$$', name(1)), 'Interpreter', 'latex', 'FontSize', titlefont)
    %legend(sprintf('$$W(t)_%s$$', name(1)), 'Location', 'southwest', 'Interpreter', 'latex', 'FontSize', labelfont)

    ax = gca; ax.YAxis(1).Color = [0, 0, 0 ]; ax.YAxis(2).Color = colors(i);
    xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
end

print(f_noSS, '../Figures/LearningWithoutScaling.png', '-dpng', '-r300')
close(f_noSS)

%% Adding synaptic scaling:
STDP = struct('window', @Kempter1999Window, 'Kupdate', @(K, W) K + W, 'w_i', 1.0e-5, 'w_o', - 1.0475*1.0e-5);
plastopts = struct('SP', STDP, 'SS', true);

[t, thetas_full, K, Kmeans] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);
drawthetas = spikesNaN(thetas_full);
z = orderparameter(thetas_full);

figure; hold on; box on;
yyaxis left
plot(t, drawthetas, '-', 'LineWidth', 0.1, 'Color', [0, 0, 1, 1/pars.N], 'HandleVisibility', 'off')
plot(t, abs(z), '-k', 'LineWidth', 2)
ylabel('$\theta_i$','Interpreter','latex', 'FontSize', labelfont)
ylim([-pi, pi]);

yyaxis right
plot(t, Kmeans, 'LineWidth', 2, 'Color', "#0072BD")
ylabel('$\langle k \rangle$','Interpreter','latex', 'FontSize', labelfont)

ax = gca; ax.YAxis(1).Color = [0, 0, 1]; ax.YAxis(2).Color = "#0072BD";
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)

%%
degrees_in = sum(K,2);
degrees_out = sum(K,1);

figure; hold on; box on;
histogram(degrees_in, 50);
histogram(degrees_out, 50);
legend("$$k_i^{\rm in}$$", "$$k_j^{\rm out}$$", 'Interpreter','latex', 'FontSize', labelfont)
