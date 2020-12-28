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
h = 0.01; tnow = h; tend = 10;

pars.N = 10;
pars.a_n = 0.666666666666666666667;
seed = 2; rng(seed);
IC = linspace(0, 2*pi - (2*pi)/(pars.N),pars.N)';
pars.e = 0; %randcauchy(seed, pars.eta0, pars.delta, pars.N);

KMAX = 10; etaMAX = 10;

%% 1. Kempter window, no IP
if true
STDP = struct('window', @Kempter1999Window, 'Kupdate', @(K, W) K + W, 'w_i', 1.0e-5, 'w_o', - 1.0475*1.0e-5);
plastopts = struct('SP', STDP, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);
drawthetas = spikesNaN(thetas_full);
z = orderparameter(thetas_full);

STDPfigure(pars, plastopts, t, drawthetas, K, Kmeans, titlefont, labelfont, 'STDPKempter', '#0072BD');
disp('Made STDPKempter figure')
end

%% 2. Song window, no IP
if true
STDP = struct('window', @Song2017Window, 'Kupdate', @(K, W) K + KMAX*W);
plastopts = struct('SP', STDP, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);
drawthetas = spikesNaN(thetas_full);
z = orderparameter(thetas_full);

STDPfigure(pars, plastopts, t, drawthetas, K, Kmeans, titlefont, labelfont, 'STDPSong', '#D95319');
disp('Made STDPSong figure')
end

%% 3. ChrollCannon window, no IP
if true
STDP = struct('window', @ChrolCannon2012Window, 'Kupdate', @(K, W) K + KMAX*W);
plastopts = struct('SP', STDP, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);
drawthetas = spikesNaN(thetas_full);
z = orderparameter(thetas_full);

STDPfigure(pars, plastopts, t, drawthetas, K, Kmeans, titlefont, labelfont, 'STDPChrolCannon', '#77AC30');
disp('Made STDPChrolCannon figure')
end

%% 4. Kempter window, with IP
if true
STDP = struct('window', @Kempter1999Window, 'Kupdate', @(K, W) K + W, 'w_i', 1.0e-5, 'w_o', - 1.0475*1.0e-5);
plastopts = struct('SP', STDP, 'IP', true, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);
drawthetas = spikesNaN(thetas_full);
z = orderparameter(thetas_full);

STDPandIPfigure(pars, plastopts, t, drawthetas, K, Kmeans, titlefont, labelfont, 'STDPandIPKempter', '#0072BD');
disp('Made STDPandIPKempter figure')
end

%% 5. Song window, with IP
if true
STDP = struct('window', @Song2017Window, 'Kupdate', @(K, W) K + KMAX*W);
plastopts = struct('SP', STDP, 'IP', true, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);
drawthetas = spikesNaN(thetas_full);
z = orderparameter(thetas_full);

STDPandIPfigure(pars, plastopts, t, drawthetas, K, Kmeans, titlefont, labelfont, 'STDPandIPSong', '#D95319');
disp('Made STDPandIPSong figure')
end

%% 6. ChrollCannon window, with IP
if true
STDP = struct('window', @ChrolCannon2012Window, 'Kupdate', @(K, W) K + KMAX*W);
plastopts = struct('SP', STDP, 'IP', true, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);
drawthetas = spikesNaN(thetas_full);
z = orderparameter(thetas_full);

STDPandIPfigure(pars, plastopts, t, drawthetas, K, Kmeans, titlefont, labelfont, 'STDPandIPChrolCannon', '#77AC30');
disp('Made STDPandIPChrolCannon figure')
end

%% Plotting the results

function fighandle = STDPfigure(pars, plastopts, t, thetas, K, Kmeans, titlefont, labelfont, figname, color)
fighandle = figure('Renderer', 'painters', 'Position', [0, 2000, 300, 1400]); hold on; box on;

subplot(5,1,1); hold on; axis square; box on;
title('Simulation results', 'FontSize', titlefont)
yyaxis left
rasterplot(t, thetas, labelfont, 0.1);

yyaxis right
plot(t, Kmeans(2,:), 'LineWidth', 2, 'Color', color)
ax = gca; ax.YAxis(1).Color = [0, 0, 1]; ax.YAxis(2).Color = "#0072BD";
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
ylabel('$\langle k \rangle$','Interpreter','latex', 'FontSize', labelfont)
ax = gca; ax.YAxis(1).Color = [0, 0, 0];

subplot(5,1,2); hold on; axis square; box on;
title('{ \boldmath $K_{ij} $}', 'Interpreter', 'latex', 'FontSize', titlefont)
xlim([0, pars.N]); ylim([0, pars.N]);
xlabel('Presynaptic neuron j', 'FontSize', labelfont)
ylabel('Postynaptic neuron i', 'FontSize', labelfont)
imagesc(K); colormap(gray);
set(gca,'YDir','reverse');

subplot(5,1,3); hold on; axis square; axis on; box on;
histogram(K, 'Normalization', 'pdf', 'FaceColor', color, 'EdgeColor', color)
title('Synaptic strength', 'FontSize', titlefont)
xlabel('$K_{ij}$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Density', 'FontSize', labelfont)
% set(gca,'XTick', [-plastopts.KMAX, 0, plastopts.KMAX], 'XTickLabel',{'-$$K{\rm max}$$','0','$$K{\rm max}$$'}, 'TickLabelInterpreter', 'latex')

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


function fighandle = STDPandIPfigure(pars, plastopts, t, thetas, K, Kmeans, titlefont, labelfont, figname, color)
fighandle = figure('Renderer', 'painters', 'Position', [0, 2000, 300, 1400]); hold on; box on;

subplot(6,1,1); hold on; axis square; box on;
title('Simulation results', 'FontSize', titlefont)
yyaxis left
rasterplot(t, thetas, labelfont, 0.1);

yyaxis right
plot(t, Kmeans(2,:), 'LineWidth', 2, 'Color', color)
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
histogram(K, 'Normalization', 'pdf', 'FaceColor', color, 'EdgeColor', color)
title('Synaptic strength', 'FontSize', titlefont)
xlabel('$K_{ij}$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Density', 'FontSize', labelfont)
% set(gca,'XTick', [-plastopts.KMAX, 0, plastopts.KMAX], 'XTickLabel',{'-$$K{\rm max}$$','0','$$K{\rm max}$$'}, 'TickLabelInterpreter', 'latex', 'FontSize', labelfont)

subplot(6,1,4); hold on; axis square; axis on; box on;
histogram(pars.e, 'Normalization', 'pdf', 'FaceColor', color, 'EdgeColor', color)
title('Excitability', 'FontSize', titlefont)
xlabel('$\eta_{i}$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Density', 'FontSize', labelfont)
% set(gca,'XTick', [-plastopts.etaMAX, 0, plastopts.etaMAX], 'XTickLabel',{'-$$\eta{\rm max}$$','0','$$\eta{\rm max}$$'}, 'TickLabelInterpreter', 'latex', 'FontSize', labelfont)

subplot(6,1,5); hold on; axis square; box on;
title('Degree distributions', 'FontSize', titlefont)
degrees_i = sum(abs(K),2);
degrees_o = sum(abs(K),1);
histogram(degrees_i, 'Normalization', 'pdf')
histogram(degrees_o, 'Normalization', 'pdf')
legend('\boldmath$k^{\rm in}$', '\boldmath$k^{\rm out}$', 'Interpreter', 'latex', 'FontSize', labelfont)
xlabel('$k$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Density', 'FontSize', labelfont)

subplot(6,1,6); hold on; axis square; box on;
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
