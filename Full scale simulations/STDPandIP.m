clear all; close all; clc;
% In this script we will validate some of the claims that have been made in
% Song2017, based on the theta neuron model.

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 15;
labelfont = 15;
export = true;

%% Make a GPU init handle:
if gpuDeviceCount > 0
    d = gpuDevice(gpuDeviceCount-1);
    disp(d)
end
initarray = make_GPUhandle();

%% Theta model parameters:
h = 0.01; tnow = h; tend = 100;

pars.N = 50;
pars.a_n = 0.666666666666666666667;
seed = 2; rng(seed);
IC = linspace(0, 2*pi - (2*pi)/(pars.N),pars.N)';
pars.e = 0; %randcauchy(seed, pars.eta0, pars.delta, pars.N);

%% 1. Kempter window, no IP
KMAX = 10; etaMAX = 10;
STDP = struct('window', @Kempter1999Window, 'Kupdate', @(K, W) K + W, 'w_i', 1.0e-5, 'w_o', - 1.0475*1.0e-5);
STDP = struct('window', @Song2017Window, 'Kupdate', @(K, W) K + KMAX.*W);
IP = struct(
plastopts = struct('SP', STDP, 'IP', IP, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);
drawthetas = spikesNaN(thetas_full);
z = orderparameter(thetas_full);

STDPfigure(pars, plastopts, t, drawthetas, K, Kmeans, titlefont, labelfont, 'SongSTDP')

%% Plotting the results

function fighandle = STDPfigure(pars, plastopts, t, thetas, K, Kmeans, titlefont, labelfont, figname)
fighandle = figure('Renderer', 'painters', 'Position', [0, 2000, 300, 1400]); hold on; box on;

subplot(5,1,1); hold on; axis square; box on;
title('Simulation results', 'FontSize', titlefont)
yyaxis left
rasterplot(t, thetas, labelfont, 0.1);

yyaxis right
plot(t, Kmeans(2,:), 'LineWidth', 2, 'Color', "#0072BD")
ax = gca; ax.YAxis(1).Color = [0, 0, 1]; ax.YAxis(2).Color = "#0072BD";
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
ylabel('$\langle k \rangle$','Interpreter','latex', 'FontSize', labelfont)
ax = gca; ax.YAxis(1).Color = [0, 0, 0];

subplot(5,1,2); hold on; axis square; box on;
title('{ \boldmath $K_{ij} $}', 'Interpreter', 'latex', 'FontSize', titlefont)
xlim([0, pars.N]); ylim([0, pars.N]);
xlabel('Presynaptic neuron j', 'FontSize', labelfont)
ylabel('Postynaptic neuron i', 'FontSize', labelfont)
imagesc(K/plastopts.KMAX); colormap(gray);
set(gca,'YDir','reverse');

subplot(5,1,3); hold on; axis square; axis on; box on;
histogram(K/plastopts.KMAX, 'Normalization', 'pdf')
title('Connectivity strength', 'FontSize', titlefont)
xlabel('$K_{ij}/K^{\rm max}$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Density', 'FontSize', labelfont)

subplot(5,1,4); hold on; axis square; box on;
title('Degree distributions', 'FontSize', titlefont)
degrees_i = sum(abs(K),2);
degrees_o = sum(abs(K),1);
histogram(degrees_i, 'Normalization', 'pdf')
histogram(degrees_o, 'Normalization', 'pdf')

legend('\boldmath$k^{\rm in}$', '\boldmath$k^{\rm out}$', 'Interpreter', 'latex', 'FontSize', labelfont)
xlabel('$k$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Density', 'FontSize', labelfont)

subplot(5,1,5); hold on; axis square; box on;
title('{ \boldmath $k^{\rm in} \leftrightarrow k^{\rm out}$ }', 'Interpreter', 'latex', 'FontSize', titlefont)
scatter(degrees_i, degrees_o, '.k')
xlabel('\boldmath$k^{\rm in}$','Interpreter','latex', 'FontSize', labelfont)
ylabel('\boldmath$k^{\rm out}$','Interpreter','latex', 'FontSize', labelfont)

MP=get(0,'MonitorPositions');
if size(MP,1)>1
    pos=get(fighandle,'Position');
    pause(0.01); % this seems sometimes necessary on a Mac
    set(fighandle,'Position',[pos(1,2)+MP(2,1:2) pos(3:4)]);
end

print(fighandle, ['../Figures/Learning/', figname, '.png'], '-dpng', '-r400')

end


function fighandle = STDPandIPfigure(pars, plastopts, t, thetas, K, Kmeans, titlefont, labelfont)
fighandle = figure('Renderer', 'painters', 'Position', [0, 2000, 300, 1400]); hold on; box on;

subplot(6,1,1); hold on; axis square; box on;
title('Simulation results', 'FontSize', titlefont)
yyaxis left
rasterplot(t, thetas, labelfont, 0.1);

yyaxis right
plot(t, Kmeans(2,:), 'LineWidth', 2, 'Color', "#0072BD")
ax = gca; ax.YAxis(1).Color = [0, 0, 1]; ax.YAxis(2).Color = "#0072BD";
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
ylabel('$\langle k \rangle$','Interpreter','latex', 'FontSize', labelfont)
ax = gca; ax.YAxis(1).Color = [0, 0, 0];

subplot(6,1,2); hold on; axis square; box on;
title('{ \boldmath $K_{ij} $}', 'Interpreter', 'latex', 'FontSize', titlefont)
xlim([0, pars.N]); ylim([0, pars.N]);
xlabel('Presynaptic neuron j', 'FontSize', labelfont)
ylabel('Postynaptic neuron i', 'FontSize', labelfont)
imagesc(K/plastopts.KMAX); colormap(gray);
set(gca,'YDir','reverse');

subplot(6,1,3); hold on; axis square; axis on; box on;
histogram(K/plastopts.KMAX, 'Normalization', 'pdf')
title('Connectivity strength', 'FontSize', titlefont)
xlabel('$K_{ij}/K^{\rm max}$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Density', 'FontSize', labelfont)

subplot(6,1,4); hold on; axis square; box on;
title('Degree distributions', 'FontSize', titlefont)
degrees_i = sum(abs(K),2);
degrees_o = sum(abs(K),1);
histogram(degrees_i, 'Normalization', 'pdf')
histogram(degrees_o, 'Normalization', 'pdf')

legend('\boldmath$k^{\rm in}$', '\boldmath$k^{\rm out}$', 'Interpreter', 'latex', 'FontSize', labelfont)
xlabel('$k$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Density', 'FontSize', labelfont)

subplot(6,1,5); hold on; axis square; box on;
title('{ \boldmath $k^{\rm in} \leftrightarrow k^{\rm out}$ }', 'Interpreter', 'latex', 'FontSize', titlefont)
scatter(degrees_i, degrees_o, '.k')
xlabel('\boldmath$k^{\rm in}$','Interpreter','latex', 'FontSize', labelfont)
ylabel('\boldmath$k^{\rm out}$','Interpreter','latex', 'FontSize', labelfont)

MP=get(0,'MonitorPositions');
if size(MP,1)>1
    pos=get(fighandle,'Position');
    pause(0.01); % this seems sometimes necessary on a Mac
    set(fighandle,'Position',[pos(1,2)+MP(2,1:2) pos(3:4)]);
end

% print(fighandle, '../Figures/STDPbeforeIP.png', '-dpng', '-r400')

end
