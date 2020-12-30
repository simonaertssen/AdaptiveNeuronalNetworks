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

pars.N = 10;
pars.a_n = 0.666666666666666666667;
seed = 2; rng(seed);
IC = linspace(0, 2*pi - (2*pi)/(pars.N),pars.N)';
pars.e = zeros(pars.N, 1); %randcauchy(seed, pars.eta0, pars.delta, pars.N);

KMAX = 10; etaMAX = 10;

%% 1. Kempter window, no IP
if true
STDP = struct('window', @Kempter1999Window, 'Kupdate', @(K, W) K + W, 'w_i', 1.0e-5, 'w_o', - 1.0475*1.0e-5);
plastopts = struct('SP', STDP, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);
% drawthetas = spikesNaN(thetas_full);

STDPfigure(pars, plastopts, t, thetas_full, K, Kmeans, titlefont, labelfont, 'STDPKempter', '#0072BD', export);
end

%% 2. Song window, no IP
if true
STDP = struct('window', @Song2017Window, 'Kupdate', @(K, W) K + KMAX*W);
plastopts = struct('SP', STDP, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);
drawthetas = spikesNaN(thetas_full);
z = orderparameter(thetas_full);

STDPfigure(pars, plastopts, t, thetas_full, K, Kmeans, titlefont, labelfont, 'STDPSong', '#D95319', export);
end

%% 3. ChrollCannon window, no IP
if true
STDP = struct('window', @ChrolCannon2012Window, 'Kupdate', @(K, W) K + KMAX*W);
plastopts = struct('SP', STDP, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);

STDPfigure(pars, plastopts, t, thetas_full, K, Kmeans, titlefont, labelfont, 'STDPChrolCannon', '#77AC30', export);
end

%% 4. Kempter window, with IP
if true
STDP = struct('window', @Kempter1999Window, 'Kupdate', @(K, W) K + W, 'w_i', 1.0e-5, 'w_o', - 1.0475*1.0e-5);
plastopts = struct('SP', STDP, 'IP', true, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);

STDPfigure(pars, plastopts, t, thetas_full, K, Kmeans, titlefont, labelfont, 'STDPandIPKempter', '#0072BD', export);
end

%% 5. Song window, with IP
if true
STDP = struct('window', @Song2017Window, 'Kupdate', @(K, W) K + KMAX*W);
plastopts = struct('SP', STDP, 'IP', true, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);

STDPfigure(pars, plastopts, t, thetas_full, K, Kmeans, titlefont, labelfont, 'STDPandIPSong', '#D95319', export);
end

%% 6. ChrollCannon window, with IP
if true
STDP = struct('window', @ChrolCannon2012Window, 'Kupdate', @(K, W) K + KMAX*W);
plastopts = struct('SP', STDP, 'IP', true, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);

STDPfigure(pars, plastopts, t, thetas_full, K, Kmeans, titlefont, labelfont, 'STDPandIPChrolCannon', '#77AC30', export);
end

%% Plotting the results
clc; close all;

% STDPfigure(pars, plastopts, t, thetas_full, K, Kmeans, titlefont, labelfont, 'test', '#77AC30', true);

% STDPfigure(pars, plastopts, t, thetas_full, K, Kmeans, titlefont, labelfont, 'test', '#D95319', true);

function fighandle = STDPfigure(pars, plastopts, t, thetas, K, Kmeans, titlefont, labelfont, figname, color, export)
numfigs = 5;
e = 0;
if isfield(plastopts, 'IP')
    numfigs = 6;
    e = 1;
end

fighandle = figure('Renderer', 'painters', 'Position', [0, 2000, 300, 1400]); hold on; box on;

sbplt(1) = subplot(numfigs,1,1); hold on; axis square; box on;
title('Simulation results', 'FontSize', titlefont, 'FontWeight', 'normal')
yyaxis left; ylim([0, 1]); xlim([t(1), t(end)]); 
z = orderparameter(thetas);
% windowSize = 500;
% b = (1/windowSize)*ones(1,windowSize);
b = normpdf(-3:0.005:3, 0, 1);
b = b/sum(b);
a = 1;
zfilt = filter(b,a,z);
ylabel('$Z(t)$','Interpreter','latex', 'FontSize', labelfont)
plot(t, abs(zfilt), '-k', 'LineWidth', 2)
% rasterplot(t, thetas, labelfont, 0.1);
yyaxis right
plot(t, Kmeans(2,:), 'LineWidth', 2, 'Color', color)
ax = gca; ax.YAxis(1).Color = [0, 0, 1]; ax.YAxis(2).Color = color;
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
ylabel('$\langle k \rangle$','Interpreter','latex', 'FontSize', labelfont)
ax = gca; ax.YAxis(1).Color = [0, 0, 0];
pos = sbplt(1).Position;
w = 0.09; y = 0.08/16;
sbplt(1).Position = [pos(1), pos(2) + w/2 + 1.5*e*y, pos(3), pos(4)];

sbplt(2) = subplot(numfigs,1,2); hold on; axis square; box on;
title('{ \boldmath $K_{ij} $}', 'Interpreter', 'latex', 'FontSize', titlefont, 'FontWeight', 'normal')
xlim([0, pars.N]); ylim([0, pars.N]);
xlabel('Presynaptic neuron j', 'FontSize', labelfont)
ylabel('Postynaptic neuron i', 'FontSize', labelfont)
imagesc(K); colormap(gray);
set(gca,'YDir','reverse');
colorbar('Location', 'southoutside')
pos = sbplt(2).Position;
sbplt(2).Position = [pos(1) - w/2, pos(2) + y + w/2, pos(3) + w, pos(4)];

sbplt(3) = subplot(numfigs,1,3); hold on; axis square; axis on; box on;
histogram(K, 'Normalization', 'pdf', 'FaceColor', color, 'FaceAlpha', 1)
title('Synaptic strength', 'FontSize', titlefont, 'FontWeight', 'normal')
xlabel('$K_{ij}$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Density', 'FontSize', labelfont)
% set(gca,'XTick', [-plastopts.KMAX, 0, plastopts.KMAX], 'XTickLabel',{'-$$K{\rm max}$$','0','$$K{\rm max}$$'}, 'TickLabelInterpreter', 'latex', 'FontSize', labelfont)
pos = sbplt(3).Position;
sbplt(3).Position = [pos(1), pos(2) + w/4 - 1.5*e*y, pos(3), pos(4)];

sbplt(4) = subplot(numfigs,1,4); hold on; axis square; box on;
title('Degree distributions', 'FontSize', titlefont, 'FontWeight', 'normal')
degrees_i = sum(abs(K),2)';
degrees_o = sum(abs(K),1);
mindegree = min([degrees_i, degrees_o]);
maxdegree = max([degrees_i, degrees_o]);
dist = maxdegree - mindegree; perc = 0.05*dist;
numbins = round(sqrt(pars.N));
xlim([mindegree - perc, maxdegree + perc]);
histogram(degrees_i, numbins, 'Normalization', 'pdf', 'FaceColor', '#000000', 'FaceAlpha', 0.3)
histogram(degrees_o, numbins, 'Normalization', 'pdf', 'FaceColor', '#d8d778', 'FaceAlpha', 0.7)
% leg = legend('\boldmath$k^{\rm in}$', '\boldmath$k^{\rm out}$', 'Interpreter', 'latex', 'FontSize', labelfont-6, 'Location', 'northwest', 'Orientation','vertical');
minx = sbplt(4).XLim(1); maxy = sbplt(4).YLim(end);
t1 = text(1.1*minx, 0.9*maxy, '\boldmath$k^{\rm in}$', 'Interpreter', 'latex', 'FontSize', labelfont-2, 'Color', '#000000');
t2 = text(t1.Position(1) + 1.5*t1.Extent(3), 0.9*maxy, '\boldmath$k^{\rm out}$', 'Interpreter', 'latex', 'FontSize', labelfont-2, 'Color', '#d8d778');
xlabel('$k$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Density', 'FontSize', labelfont)
pos = sbplt(4).Position;
h = 0.02;
sbplt(4).Position = [pos(1), pos(2) + h - 2*e*y, pos(3), pos(4)];

sbplt(5) = subplot(numfigs,1,5); hold on; axis square; box on;
xlim([mindegree - perc, maxdegree + perc]);
ylim([mindegree - perc, maxdegree + perc]);
% title('{ \boldmath $k^{\rm in} \leftrightarrow k^{\rm out}$ }', 'Interpreter', 'latex', 'FontSize', titlefont)
scatter(degrees_i, degrees_o, 15, 'MarkerEdgeColor', '#000000', 'MarkerFaceColor', color)
xlabel('\boldmath$k^{\rm in}$','Interpreter','latex', 'FontSize', labelfont)
ylabel('\boldmath$k^{\rm out}$','Interpreter','latex', 'FontSize', labelfont)
pos = sbplt(5).Position;
sbplt(5).Position = [pos(1), pos(2) + 2*h - 4*e*y, pos(3), pos(4)];

if isfield(plastopts, 'IP')
sbplt(6) = subplot(numfigs,1,6); hold on; axis square; axis on; box on;
histogram(pars.e, 'Normalization', 'pdf', 'FaceColor', color, 'FaceAlpha', 1)
title('Excitability', 'FontSize', titlefont, 'FontWeight', 'normal')
xlabel('$\eta_{i}$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Density', 'FontSize', labelfont)
% set(gca,'XTick', [-plastopts.etaMAX, 0, plastopts.etaMAX], 'XTickLabel',{'-$$\eta{\rm max}$$','0','$$\eta{\rm max}$$'}, 'TickLabelInterpreter', 'latex', 'FontSize', labelfont)
pos = sbplt(6).Position;
sbplt(6).Position = [pos(1), pos(2) + 2.5*h - 7*e*y, pos(3), pos(4)];
end


MP=get(0,'MonitorPositions');
if size(MP,1)>1
    pos=get(fighandle,'Position');
    pause(0.01); % this seems sometimes necessary on a Mac
    set(fighandle,'Position',[pos(1,2)+MP(2,1:2) pos(3:4)]);
end

set(findall(gcf,'-property','FontName'),'FontName','Avenir')

if export
exportgraphics(fighandle,['../Figures/Learning/', figname, '.pdf'], 'ContentType','vector')
end
disp(['Made ', figname, ' figure'])
end
