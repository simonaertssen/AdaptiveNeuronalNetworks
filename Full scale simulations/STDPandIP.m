clear all; close all; clc;
% In this script we will validate some of the claims that have been made in
% Song2017, based on the theta neuron model.

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 18;
labelfont = 15;
export = true;

%% Make a GPU init handle:
if gpuDeviceCount > 0
    d = gpuDevice(gpuDeviceCount-1);
    disp(d)
end
initarray = make_GPUhandle();

%% Theta model parameters:
h = 0.01; tnow = 0; tend = 4000;

pars.N = 100;
pars.a_n = 0.666666666666666666667;
seed = 2; rng(seed);
IC = linspace(0, 2*pi - (2*pi)/(pars.N),pars.N)';
pars.e = zeros(pars.N, 1); %randcauchy(seed, pars.eta0, pars.delta, pars.N);

KMAX = 100; etaMAX = 100;

%% Only STDP:
if true
fighandle = figure('Renderer', 'painters', 'Position', [0, 2000, 800, 1400]); 

%% 1. Kempter window, no IP
STDP = struct('window', @Kempter1999Window, 'Kupdate', @(K, W) K + W, 'w_i', 1.0e-3, 'w_o', - 1.0475*1.0e-3);
plastopts = struct('SP', STDP);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);
% drawthetas = spikesNaN(thetas_full);

STDPfigure(0, pars, plastopts, t, thetas_full, K, Kmeans, titlefont, labelfont, 'STDPKempter', '#298A3E', export);

%% 2. Song window, no IP
STDP = struct('window', @Song2012Window, 'Kupdate', @(K, W) K + KMAX*W);
plastopts = struct('SP', STDP, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);

STDPfigure(1, pars, plastopts, t, thetas_full, K, Kmeans, titlefont, labelfont, 'STDPSong', '#D95319', export);

%% 3. ChrollCannon window, no IP
STDP = struct('window', @ChrolCannon2012Window, 'Kupdate', @(K, W) K + KMAX*W);
plastopts = struct('SP', STDP, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);

STDPfigure(2, pars, plastopts, t, thetas_full, K, Kmeans, titlefont, labelfont, 'STDPChrolCannon', '#EDB120', export);
%% Move figure to big screen:
moveFigToBigScreen(fighandle);
set(findall(gcf,'-property','FontName'),'FontName','Avenir')

%% Export figure:
if export
exportgraphics(fighandle,'../Figures/Learning/STDP.pdf', 'ContentType','vector')
end
close(fighandle)
disp('Made STDP figure')

end

%% STDP and IP figure:
if true
fighandle = figure('Renderer', 'painters', 'Position', [0, 2000, 800, 1400]); 

%% 4. Kempter window, with IP
STDP = struct('window', @Kempter1999Window, 'Kupdate', @(K, W) K + W, 'w_i', 1.0e-3, 'w_o', - 1.0475*1.0e-3);
plastopts = struct('SP', STDP, 'IP', true, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);

STDPfigure(0, pars, plastopts, t, thetas_full, K, Kmeans, titlefont, labelfont, 'STDPandIPKempter', '#298A3E', export);

%% 5. Song window, with IP
STDP = struct('window', @Song2012Window, 'Kupdate', @(K, W) K + KMAX*W);
plastopts = struct('SP', STDP, 'IP', true, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);

STDPfigure(1, pars, plastopts, t, thetas_full, K, Kmeans, titlefont, labelfont, 'STDPandIPSong', '#D95319', export);

%% 6. ChrollCannon window, with IP
STDP = struct('window', @ChrolCannon2012Window, 'Kupdate', @(K, W) K + KMAX*W);
plastopts = struct('SP', STDP, 'IP', true, 'KMAX', KMAX, 'etaMAX', etaMAX);

[t, thetas_full, K, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);

STDPfigure(2, pars, plastopts, t, thetas_full, K, Kmeans, titlefont, labelfont, 'STDPandIPChrolCannon', '#EDB120', export);

%% Move figure to big screen:
set(findall(gcf,'-property','FontName'),'FontName','Avenir')
moveFigToBigScreen(fighandle);
end

if export
exportgraphics(fighandle,'../Figures/Learning/STDPandIP.pdf', 'ContentType','vector')
end
disp('Made STDP and IP figure')
close(fighandle)

%% Plotting the results

% STDPfigure(0, pars, plastopts, t, thetas_full, K, Kmeans, titlefont, labelfont, 'STDPKempter', '#0072BD', export);

function fighandle = STDPfigure(idx, pars, plastopts, t, thetas, K, Kmeans, titlefont, labelfont, figname, color, export)
numfigs = 5;
e = 0;
if isfield(plastopts, 'IP')
    numfigs = 6;
    e = 1;
end

sbplt(1) = subplot(numfigs,3,1+idx); hold on; axis square; box on;
title('Simulation results', 'FontSize', titlefont, 'FontWeight', 'normal')
yyaxis left; ylim([0, 1]); xlim([t(1), t(end)]); 
z = orderparameter(thetas);
b = normpdf(-3:0.005:3, 0, 1); b = b/sum(b); a = 1; zfilt = filter(b,a,z);
if idx == 0; ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont); end
plot(t, abs(zfilt), '-k', 'LineWidth', 2)
xticks(linspace(0, t(end), 5))

yyaxis right
plot(t, Kmeans(2,:), 'LineWidth', 2, 'Color', color)
ax = gca; ax.YAxis(1).Color = [0, 0, 0]; ax.YAxis(2).Color = color;
if abs(ax.YTick(2) - ax.YTick(1)) < 1
    ax.YTick = sort([Kmeans(2,1), Kmeans(2,end)]);
    ax.YTickLabel = string([floor(Kmeans(2,1)), ceil(Kmeans(2,end))]);
end
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
if idx == 2; ylabel('$\langle k \rangle$','Interpreter','latex', 'FontSize', labelfont); end
pos = sbplt(1).Position;
w = 0.09; y = 0.08/16;
sbplt(1).Position = [pos(1), pos(2) + w/1.7 + 1.5*e*y, pos(3), pos(4)];

sbplt(2) = subplot(numfigs,3,4+idx); hold on; axis square; box on;
title('{ \boldmath $K_{ij} $}', 'Interpreter', 'latex', 'FontSize', titlefont, 'FontWeight', 'normal')
xlim([0, pars.N]); ylim([0, pars.N]);
xlabel('Presynaptic neuron j', 'FontSize', labelfont)
if idx == 0; ylabel('Postynaptic neuron i', 'FontSize', labelfont); end
im = imagesc(K); colormap(gray);
xdata = im.XData; ydata = im.YData;
im.XData = xdata - 0.5; im.YData = ydata - 0.5;
set(gca,'YDir','reverse');
colorbar('Location', 'southoutside')
pos = sbplt(2).Position;
sbplt(2).Position = [pos(1) - w/2, pos(2) + y + w/2 + 0.5*e*y, pos(3) + w, pos(4)];

sbplt(3) = subplot(numfigs,3,7+idx); hold on; axis square; axis on; box on;
histogram(K, 'Normalization', 'pdf', 'FaceColor', color, 'FaceAlpha', 1)
title('Synaptic strength', 'FontSize', titlefont, 'FontWeight', 'normal')
xlabel('$K_{ij}$','Interpreter','latex', 'FontSize', labelfont)
if idx == 0; ylabel('Density', 'FontSize', labelfont); end
pos = sbplt(3).Position;
sbplt(3).Position = [pos(1), pos(2) + w/4 - 1.5*e*y, pos(3), pos(4)];

sbplt(4) = subplot(numfigs,3,10+idx); hold on; axis square; box on;
title('Degree distributions', 'FontSize', titlefont, 'FontWeight', 'normal')
degrees_i = sum(abs(K),2)';
degrees_o = sum(abs(K),1);
mindegree = min([degrees_i, degrees_o]);
maxdegree = max([degrees_i, degrees_o]);
dist = maxdegree - mindegree; perc = 0.05*dist;
numbins = round(sqrt(pars.N));
xlim([mindegree - perc, maxdegree + perc]);
histogram(degrees_i, linspace(mindegree, maxdegree, numbins), 'Normalization', 'pdf', 'FaceColor', '#000000', 'FaceAlpha', 0.3)
histogram(degrees_o, linspace(mindegree, maxdegree, numbins), 'Normalization', 'pdf', 'FaceColor', '#d8d778', 'FaceAlpha', 0.7)
minx = sbplt(4).XLim(1); maxy = sbplt(4).YLim(end); dist = sbplt(4).XLim(end) - minx;
ylim([0, 1.06*maxy]);
t1 = text(minx + 0.02*dist, 1.04*maxy, '\boldmath$k^{\rm in}$', 'Interpreter', 'latex', 'FontSize', labelfont-2, 'Color', '#000000', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
t2 = text(t1.Position(1) + 1.5*t1.Extent(3), 1.04*maxy, '\boldmath$k^{\rm out}$', 'Interpreter', 'latex', 'FontSize', labelfont-2, 'Color', '#d8d778', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
xlabel('$k$','Interpreter','latex', 'FontSize', labelfont)
if idx == 0; ylabel('Density', 'FontSize', labelfont); end
pos = sbplt(4).Position;
h = 0.02;
sbplt(4).Position = [pos(1), pos(2) + h - 2*e*y, pos(3), pos(4)];

sbplt(5) = subplot(numfigs,3,13+idx); hold on; axis square; box on;
xlim([mindegree - perc, maxdegree + perc]);
ylim([mindegree - perc, maxdegree + perc]);
% title('{ \boldmath $k^{\rm in} \leftrightarrow k^{\rm out}$ }', 'Interpreter', 'latex', 'FontSize', titlefont)
scatter(degrees_i, degrees_o, 15, 'MarkerEdgeColor', '#000000', 'MarkerFaceColor', color)
xlabel('\boldmath$k^{\rm in}$','Interpreter','latex', 'FontSize', labelfont)
if idx == 0; ylabel('\boldmath$k^{\rm out}$','Interpreter','latex', 'FontSize', labelfont); end
pos = sbplt(5).Position;
sbplt(5).Position = [pos(1), pos(2) + 2*h - 4*e*y, pos(3), pos(4)];

if isfield(plastopts, 'IP')
sbplt(6) = subplot(numfigs,3,16+idx); hold on; axis square; axis on; box on;
histogram(pars.e, linspace(-plastopts.etaMAX, plastopts.etaMAX, 21), 'Normalization', 'pdf', 'FaceColor', color, 'FaceAlpha', 1)
title('Excitability', 'FontSize', titlefont, 'FontWeight', 'normal')
xlabel('$\eta_{i}$','Interpreter','latex', 'FontSize', labelfont)
if idx == 0; ylabel('Density', 'FontSize', labelfont); end
% set(gca,'XTick', [-plastopts.etaMAX, 0, plastopts.etaMAX], 'XTickLabel',{'-$$\eta{\rm max}$$','0','$$\eta{\rm max}$$'}, 'TickLabelInterpreter', 'latex', 'FontSize', labelfont)
pos = sbplt(6).Position;
sbplt(6).Position = [pos(1), pos(2) + 2*h - 5*e*y, pos(3), pos(4)];
end

end
