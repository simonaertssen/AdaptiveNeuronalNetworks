clear all; close all; clc;
% Here we will make phase space plots for the report, illustrating all the
% different global states

%%
addpath('../Functions');
addpath('../Mean Field Reductions/');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')
cm = [1,0,0; 0, 0.7410, 0.4470; 0, 0.4470, 0.6410];

titlefont = 15;
linewidth = 4;
linewidth_total = 2.5;
arrowst = 12;
labelfont = 17;

fsolveoptions = optimset('Display','off');
odeoptions = odeset('RelTol', 1.0e-8,'AbsTol', 1.0e-8);

%% Grids:
% Axis square:
unitsquare = [-1, -1, 1, 1; 1, -1, -1, 1];
th = 0:pi/200:2*pi;
unitcircle = [cos(th); sin(th)];
drawcircle = round(unitcircle);

% Main grid:
l = 1; stp = 2*l/40; interval = -l:stp:l;
[X,Y] = meshgrid(interval,interval);

[in, ~] = inpolygon(X, Y, cos(th), sin(th));
[sz, ~] = size(in); m = ceil(sz/2); l = round(m*0.25);
in(1, :) = 0; in(end, :) = 0; in(:, 1) = 0; in(:, end) = 0;
in(1, m-l:m+l) = 1; in(end, m-l:m+l) = 1; in(m-l:m+l, 1) = 1; in(m-l:m+l, end) = 1;

%% Theta neurons parameters:
pars.N = 1000;
pars.a_n = 0.666666666666666666667;
tnow = 0; tend = 100;
seed = 3; rng(seed);

%% Draw the problems with mappings, use PSR state:
f_finalcond = figure('Renderer', 'painters', 'Position', [-500 1600 1000 200]); 
Zstart = -0.2 + 1i*0.8;

% Random network:
pars.eta0 = -0.2; pars.delta = 0.1; pars.K = -2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
p = prepareOAparameters(make_randomparameters(pars, 0.3));
idx = round(linspace(1, p.Mk, 40));

OAIC = ones(p.Mk,1)*Zstart;
[~, ZOA, bs] = OA_simulatenetwork(0, tend, OAIC, p, odeoptions);

subplot(1,4,1); hold on; box on; axis square
length = norm(bs(end,end) - bs(end,1));
xlim([min(real(bs(end,:))) - length*0.05, max(real(bs(end,:))) + length*0.05]);
ylim([min(imag(bs(end,:))) - length*0.05, max(imag(bs(end,:))) + length*0.05]);
plot(real(bs(end,:)), imag(bs(end,:)), 'Color', p.colorvec)
scatter(real(bs(end,idx)), imag(bs(end,idx)), 200, p.colorvec, '.')
ZOA = ZOA(1:370);
plot(real(ZOA), imag(ZOA), 'LineWidth', linewidth_total, 'color', p.colorvec);
endline = ZOA(end-100) - ZOA(end);
endpoint = ZOA(end) + 0.008*endline/abs(endline);
plot_arrow(real(endpoint), imag(endpoint), real(ZOA(end)), imag(ZOA(end)),'linewidth', 2, ...
    'color', p.colorvec,'facecolor', p.colorvec,'edgecolor', p.colorvec, 'headwidth',0.5,'headheight',2);
h = plot(cos(th), sin(th), '-k', 'LineWidth', 2);
text(real(bs(end,1)), imag(bs(end,2)), '1', 'FontSize', titlefont-2, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
text(real(bs(end,end-2)), imag(bs(end,end)), '0', 'FontSize', titlefont-2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

subplot(1,4,2); hold on; box on;
[edges, counts] = projectToHist(flip(bs(end,:)));
xlim([0, 1]);
histogram('BinEdges',edges,'BinCounts',counts, 'Normalization', 'pdf', 'FaceColor', p.colorvec)
xlabel('Curve length', 'FontSize', labelfont)
ylabel('Density', 'FontSize', labelfont)


% Scale-free network:
pars.eta0 = -0.2; pars.delta = 0.1; pars.K = -2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
p = prepareOAparameters(make_scalefreeparameters(pars, 3));
idx = round(linspace(1, p.Mk, 25));

OAIC = ones(p.Mk,1)*Zstart;
[~, ZOA, bs] = OA_simulatenetwork(0, tend, OAIC, p, odeoptions);

sp = subplot(1,4,3); hold on; box on; axis square
xlim([min(real(bs(end,:))) - length*0.05, max(real(bs(end,:))) + length*0.05]);
ylim([min(imag(bs(end,:))) - length*0.05, max(imag(bs(end,:))) + length*0.05]);
plot(real(bs(end,:)), imag(bs(end,:)), 'Color', p.colorvec)
scatter(real(bs(end,idx)), imag(bs(end,idx)), 200, p.colorvec, '.')
ZOA = ZOA(1:350);
plot(real(ZOA), imag(ZOA), 'LineWidth', linewidth_total, 'color', p.colorvec);
endline = ZOA(end-50) - ZOA(end);
endpoint = ZOA(end) + 0.02*endline/abs(endline);
plot_arrow(real(endpoint), imag(endpoint), real(ZOA(end)), imag(ZOA(end)),'linewidth', 2, ...
    'color', p.colorvec,'facecolor', p.colorvec,'edgecolor', p.colorvec, 'headwidth',0.5,'headheight',2);
h = plot(cos(th), sin(th), '-k', 'LineWidth', 2);
text(real(bs(end,3)), imag(bs(end,5)), '1', 'FontSize', titlefont-2, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
text(real(bs(end,end-20)), imag(bs(end,end-2)), '0', 'FontSize', titlefont-2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
% pos = sp.Position
% sp.Position = [pos(1)+0.03, pos(2)+0.015, pos(3)-0.01, pos(4)-0.01]

sp = subplot(1,4,4); hold on; box on;
[edges, counts] = projectToHist(flip(bs(end,:)));
xlim([0, 1]);
histogram('BinEdges',edges,'BinCounts',counts, 'Normalization', 'pdf', 'FaceColor', p.colorvec)
xlabel('Curve length', 'FontSize', labelfont)
ylabel('Density', 'FontSize', labelfont)
% pos = sp.Position
% sp.Position = [pos(1)+0.03, pos(2)+0.015, pos(3)-0.03, pos(4)-0.03]

% Export the figure:
set(findall(gcf,'-property','FontName'),'FontName','Avenir')
exportgraphics(f_finalcond,'../Figures/Distributions/FinalConditions.pdf')


%% Functions:
function [edges, counts] = projectToHist(endb)
    Mk = numel(endb);
    nbins = round(sqrt(Mk));
    counts = zeros(nbins,1);
    % Calculate the total difference
    differences = diff(endb);
    totaldist = 0;
    for k = 1:Mk-1; totaldist = totaldist + norm(differences(k)); end
    % Collect counts per bin
    binlength = totaldist/nbins; binidx = 1; bindist = binlength;
    for k = 1:Mk
        dist = norm(endb(k) - endb(1));
        if dist > bindist
            binidx = binidx + 1;
            bindist = bindist + binlength; 
        end
        counts(binidx) = counts(binidx) + 1;
    end
    edges = linspace(0,1,nbins+1);
end
