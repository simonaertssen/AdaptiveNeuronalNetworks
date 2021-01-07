clear all; close all; clc;
% In this script we will be investigating the influence of the size of the
% in- and output weights put on in- and outgoing signals using STDP as
% formulated by Kempter.

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 18;
labelfont = 20;
export = true;

%% Make a GPU init handle:
if gpuDeviceCount > 0
    d = gpuDevice(gpuDeviceCount-1);
    disp(d)
end
initarray = make_GPUhandle();

%% Theta model parameters:
h = 0.01; tnow = h; tend = 3000;

pars.N = 100;
pars.a_n = 0.666666666666666666667;
seed = 2; rng(seed);
IC = linspace(0, 2*pi - (2*pi)/(pars.N),pars.N)';
pars.e = zeros(pars.N, 1); %randcauchy(seed, pars.eta0, pars.delta, pars.N);

KMAX = 10; etaMAX = 10;
color = '#3AC1D6';

%% Figure handle:
f_kempter = figure('Renderer', 'painters', 'Position', [50 800 1000 200]); 

weights = [1.0e-5, 1.0e-2, 1.0e-1];
for i = 1:3
    weight = weights(i);
    sbplt(i) = subplot(1,3,i); hold on; box on;
    
    STDP = struct('window', @Kempter1999Window, 'Kupdate', @(K, W) K + W, 'w_i', weight, 'w_o', - 1.0475*weight);
    plastopts = struct('SP', STDP, 'KMAX', KMAX, 'etaMAX', etaMAX);

    [t, thetas_full, ~, Kmeans, pars] = DOPRI_simulatenetwork_adaptive(tnow,tend,IC,h,pars,plastopts);
    
    title(sprintf('$$w^{\\rm in}$$ = %0.1E', weight), 'FontSize', titlefont, 'FontWeight', 'normal', 'Interpreter','latex')
    yyaxis left; ylim([0, 1]); xlim([t(1), t(end)]);
    z = orderparameter(thetas_full);
    b = normpdf(-3:0.005:3, 0, 1); b = b/sum(b); a = 1; zfilt = filter(b,a,z);
    if i == 1; ylabel('$Z(t)$','Interpreter','latex', 'FontSize', labelfont); end
    plot(t, abs(zfilt), '-k', 'LineWidth', 2)
    
    yyaxis right
    plot(t, Kmeans(2,:), 'LineWidth', 2, 'Color', color)
    ax = gca; ax.YAxis(1).Color = [0, 0, 0]; ax.YAxis(2).Color = color;
    if abs(ax.YTick(2) - ax.YTick(1)) < 1
        ax.YTick = sort([Kmeans(2,1), Kmeans(2,end)]);
        ax.YTickLabel = string([floor(Kmeans(2,1)), ceil(Kmeans(2,end))]);
    end
    xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
    if i == 3; ylabel('$\langle k \rangle$','Interpreter','latex', 'FontSize', labelfont); end
    fprintf('Figure %d\n', i)
end

set(findall(gcf,'-property','FontName'),'FontName','Avenir')

if export
exportgraphics(f_kempter,'../Figures/Learning/KempterWinWout.pdf', 'ContentType','vector')
close(f_kempter)
end
