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

fsolveoptions = optimset('Display','off');
odeoptions = odeset('RelTol', 1.0e-8,'AbsTol', 1.0e-8);

%% Grids:
% Axis square:
unitsquare = [-1, -1, 1, 1; 1, -1, -1, 1];
th = 0:pi/50:2*pi;
unitcircle = [cos(th); sin(th)];
drawcircle = round(unitcircle);

% Main grid:
l = 1; stp = 2*l/40; interval = -l:stp:l;
[X,Y] = meshgrid(interval,interval);

[in, ~] = inpolygon(X, Y, cos(th), sin(th));
[sz, ~] = size(in); m = ceil(sz/2); l = round(m*0.25);
in(1, :) = 0; in(end, :) = 0; in(:, 1) = 0; in(:, end) = 0;
in(1, m-l:m+l) = 1; in(end, m-l:m+l) = 1; in(m-l:m+l, 1) = 1; in(m-l:m+l, end) = 1;

% Draw grid
% l = 1; stp = 2*l/6; interval = -l:stp:l;
% [X1,Y1] = meshgrid(interval,interval);
% 
% % Startpoints
% [in, on] = inpolygon(X1, Y1, [0.9*cos(th) 0], [0 0.9*sin(th)]);
% [sz, ~] = size(in); m = ceil(sz/2);
% in(1, :) = 0; in(end, :) = 0; in(:, 1) = 0; in(:, end) = 0;
% in(1, m) = 1; in(end, m) = 1; in(m, 1) = 1; in(m, end) = 1;
% X1start = reshape(X1(in), 5, 5); Y1start = reshape(Y1(in), 5, 5);

%% Theta neurons parameters:
pars.N = 1000;
pars.a_n = 0.666666666666666666667;
seed = 3; rng(seed);

%% Draw the problems with mappings, use PSR state:
labelfont = 17;
pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
p = prepareOAparameters(make_scalefreeparameters(pars, 3));

f_mappings = figure('Renderer', 'painters', 'Position', [0,0,800,800]); 
% Zstart = 0.8*cos(3*pi/5) + 1i*0.8*sin(3*pi/5);
Zstart = -0.2 + 1i*0.8; tend = 1.2; col = p.colorvec;

nlines = min(30, p.Mk); linalpha = 0.15; cols = zeros(nlines,3);
idx = randi([1 p.Mk],1,nlines);

% With the simple same IC
OAIC = ones(p.Mk,1)*Zstart;
[~, ZOA, bs] = OA_simulatenetwork(0, tend, OAIC, p, odeoptions);

% Focus on the ICs:
subplot(3, 3, 1); hold on; box on; axis square;
title('$$z($$\boldmath$$k$$\unboldmath$$, 0) = Z(0)$$', 'Interpreter', 'latex','FontSize', labelfont);
xlim([real(Zstart) - 0.1, real(Zstart) + 0.1]) 
ylim([imag(Zstart) - 0.1, imag(Zstart) + 0.1]) 
drawOAvectors(X + 1i*Y, in, p, cm(2,:));

bs_plot = plot(real(bs(:, idx)), imag(bs(:, idx)), 'LineWidth', linewidth);
for i = 1:nlines; cols(i,:) = bs_plot(i).Color; bs_plot(i).Color(4) = linalpha; end
scatter(real(Zstart), imag(Zstart), 500, [0,0,0], 'x', 'LineWidth', linewidth);
scatter(real(bs(1, idx)), imag(bs(1, idx)), 200, cols, 'x', 'LineWidth', linewidth, 'MarkerEdgeAlpha', linalpha);
bsstart_scatter = scatter(2, 2, [], cols(1,:), 'LineWidth', linewidth, 'MarkerEdgeAlpha', linalpha);
bsstart_scatter.CData(1,:) = bs_plot(1).Color;
scatter(real(ZOA(1)), imag(ZOA(1)), 500, col, '+', 'LineWidth', linewidth);
Zstart_scatter = scatter(2, 2, 500, [0,0,0], 'x', 'LineWidth', 2);
ZOAstart_scatter = scatter(2, 2, 500, col, '+', 'LineWidth', 2);
ZOA_plot = plot(real(ZOA), imag(ZOA), 'LineWidth', linewidth, 'color', p.color);

phasespaceplot();
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

% Show the final state:
subplot(3, 3, 4); hold on; box on; axis square;
xlim([-0.7, -0.3]); ylim([-0.9, -0.5]) 
drawOAvectors(X + 1i*Y, in, p, cm(2,:));

bs_plot = plot(real(bs(:, idx)), imag(bs(:, idx)), 'LineWidth', linewidth);
for i = 1:nlines; cols(i,:) = bs_plot(i).Color; bs_plot(i).Color(4) = linalpha; end
scatter(real(bs(end, idx)), imag(bs(end, idx)), 125, '.k'); 
scatter(real(bs(end,   :)), imag(bs(end,   :)), 125, '.k', 'MarkerEdgeAlpha', linalpha); 
plot(real(ZOA(1:end-10)), imag(ZOA(1:end-10)), 'LineWidth', linewidth, 'color', col);
endline = ZOA(end-3) - ZOA(end);
endpoint = ZOA(end) + 0.06*endline/abs(endline);
ZOAend_arr = plot_arrow(real(endpoint), imag(endpoint), real(ZOA(end)), imag(ZOA(end)),'linewidth', 2, ...
    'color', col,'facecolor', col,'edgecolor', col, 'headwidth',0.5,'headheight',2);

phasespaceplot();
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

% Show the full space:
subplot(3, 3, 7); hold on; box on; axis square;
xlim([-1, 1]); ylim([-1, 1])  
drawOAvectors(X + 1i*Y, in, p, cm(2,:));

bs_plot = plot(real(bs(:, idx)), imag(bs(:, idx)), 'LineWidth', linewidth_total);
for i = 1:nlines; cols(i,:) = bs_plot(i).Color; bs_plot(i).Color(4) = linalpha; end
scatter(real(bs(1, idx)), imag(bs(1, idx)), 100, cols, 'x', 'LineWidth', linewidth_total, 'MarkerEdgeAlpha', linalpha);
scatter(real(bs(end, idx)), imag(bs(end, idx)), 125, '.k'); 
scatter(real(bs(end,   :)), imag(bs(end,   :)), 125, '.k', 'MarkerEdgeAlpha', linalpha); 
scatter(real(Zstart), imag(Zstart), 200, [0,0,0], 'x', 'LineWidth', linewidth_total);
scatter(real(ZOA(1)), imag(ZOA(1)), 200, col, '+', 'LineWidth', linewidth_total);
plot(real(ZOA), imag(ZOA), 'LineWidth', linewidth_total, 'color', col);
endline = ZOA(end-4) - ZOA(end);
endpoint = ZOA(end) + 0.08*endline/abs(endline);
plot_arrow(real(endpoint), imag(endpoint), real(ZOA(end)), imag(ZOA(end)),'linewidth', 2, ...
    'color', col,'facecolor', col,'edgecolor', col, 'headwidth',0.7,'headheight',3);

phasespaceplot();
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
disp("First column")


% With the proposed analytical IC
OAIC = map_Ztozoa(conj(Zstart), p);
[~, ZOA, bs] = OA_simulatenetwork(0, tend, OAIC, p, odeoptions);

% Focus on the ICs:
subplot(3, 3, 2); hold on; box on; axis square;
title('$$z($$\boldmath$$k$$\unboldmath$$, 0) = Z(0) \cdot n($$\boldmath$$k$$\unboldmath$$) / P($$\boldmath$$k$$\unboldmath$$)$$', 'Interpreter', 'latex','FontSize', labelfont);
xlim([real(Zstart) - 0.1, real(Zstart) + 0.1]) 
ylim([imag(Zstart) - 0.1, imag(Zstart) + 0.1]) 
drawOAvectors(X + 1i*Y, in, p, cm(2,:));

bs_plot = plot(real(bs(:, idx)), imag(bs(:, idx)), 'LineWidth', linewidth);
for i = 1:nlines; cols(i,:) = bs_plot(i).Color; bs_plot(i).Color(4) = linalpha; end
scatter(real(Zstart), imag(Zstart), 500, [0,0,0], 'x', 'LineWidth', linewidth);
scatter(real(bs(1, idx)), imag(bs(1, idx)), 200, cols, 'x', 'LineWidth', linewidth, 'MarkerEdgeAlpha', linalpha);
scatter(real(ZOA(1)), imag(ZOA(1)), 500, col, '+', 'LineWidth', linewidth);
plot(real(ZOA), imag(ZOA), 'LineWidth', linewidth, 'color', col);
line([0,real(2*Zstart)], [0, imag(2*Zstart)], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k');

phasespaceplot();
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

% Show the final state:
subplot(3, 3, 5); hold on; box on; axis square;
xlim([-0.7, -0.3]); ylim([-0.9, -0.5]) 
drawOAvectors(X + 1i*Y, in, p, cm(2,:));

bs_plot = plot(real(bs(:, idx)), imag(bs(:, idx)), 'LineWidth', linewidth);
for i = 1:nlines; cols(i,:) = bs_plot(i).Color; bs_plot(i).Color(4) = linalpha; end
scatter(real(bs(end, idx)), imag(bs(end, idx)), 125, '.k'); 
scatter(real(bs(end,   :)), imag(bs(end,   :)), 125, '.k', 'MarkerEdgeAlpha', linalpha); 
plot(real(ZOA(1:end-10)), imag(ZOA(1:end-10)), 'LineWidth', linewidth, 'color', col);
endline = ZOA(end-3) - ZOA(end);
endpoint = ZOA(end) + 0.06*endline/abs(endline);
plot_arrow(real(endpoint), imag(endpoint), real(ZOA(end)), imag(ZOA(end)),'linewidth', 2, ...
    'color', col,'facecolor', col,'edgecolor', col, 'headwidth',0.5,'headheight',2);

phasespaceplot();
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

% Show the full space:
subplot(3, 3, 8); hold on; box on; axis square;
xlim([-1, 1]); ylim([-1, 1])  
drawOAvectors(X + 1i*Y, in, p, cm(2,:));

bs_plot = plot(real(bs(:, idx)), imag(bs(:, idx)), 'LineWidth', linewidth_total);
for i = 1:nlines; cols(i,:) = bs_plot(i).Color; bs_plot(i).Color(4) = linalpha; end
scatter(real(bs(1, idx)), imag(bs(1, idx)), 100, cols, 'x', 'LineWidth', linewidth_total, 'MarkerEdgeAlpha', linalpha);
scatter(real(bs(end, idx)), imag(bs(end, idx)), 100, '.k'); 
scatter(real(bs(end,   :)), imag(bs(end,   :)), 100, '.k', 'MarkerEdgeAlpha', linalpha); 
scatter(real(Zstart), imag(Zstart), 200, [0,0,0], 'x', 'LineWidth', linewidth_total);
scatter(real(ZOA(1)), imag(ZOA(1)), 200, col, '+', 'LineWidth', linewidth_total);
plot(real(ZOA), imag(ZOA), 'LineWidth', linewidth_total, 'color', col);
endline = ZOA(end-4) - ZOA(end);
endpoint = ZOA(end) + 0.08*endline/abs(endline);
plot_arrow(real(endpoint), imag(endpoint), real(ZOA(end)), imag(ZOA(end)),'linewidth', 2, ...
    'color', col,'facecolor', col,'edgecolor', col, 'headwidth',0.7,'headheight',3);
line([0,real(2*Zstart)], [0, imag(2*Zstart)], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k');

phasespaceplot();
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
disp("Second column")



% With the converged IC
[OAIC, zoaIC] = map_Ztozoa_better(Zstart,p,'fmincon');
[~, ZOA, bs] = OA_simulatenetwork(0, tend, OAIC, p, odeoptions);

% Focus on the ICs:
subplot(3, 3, 3); hold on; box on; axis square;
title('$$\hat{z}($$\boldmath$$k$$\unboldmath$$, 0) \rightarrow z($$\boldmath$$k$$\unboldmath$$, 0)$$', 'Interpreter', 'latex','FontSize', labelfont);
xlim([real(Zstart) - 0.1, real(Zstart) + 0.1]) 
ylim([imag(Zstart) - 0.1, imag(Zstart) + 0.1]) 
drawOAvectors(X + 1i*Y, in, p, cm(2,:));

bs_plot = plot(real(bs(:, idx)), imag(bs(:, idx)), 'LineWidth', linewidth);
for i = 1:nlines; cols(i,:) = bs_plot(i).Color; bs_plot(i).Color(4) = linalpha; end
endidx = randi([1 p.Mk],1,round(p.Mk/3));
scatter(real(zoaIC(endidx)), imag(zoaIC(endidx)), 25, [0,0,0], 'o');
zoa_est = scatter(2, 2, 100, [0,0,0], 'o');
scatter(real(Zstart), imag(Zstart), 500, [0,0,0], 'x', 'LineWidth', linewidth);
scatter(real(bs(1, idx)), imag(bs(1, idx)), 200, cols, 'x', 'LineWidth', linewidth, 'MarkerEdgeAlpha', linalpha);
scatter(real(ZOA(1)), imag(ZOA(1)), 500, col, '+', 'LineWidth', linewidth);
plot(real(ZOA), imag(ZOA), 'LineWidth', linewidth, 'color', col);

phasespaceplot();
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

% Show the final state:
subplot(3, 3, 6); hold on; box on; axis square;
xlim([-0.7, -0.3]); ylim([-0.9, -0.5]) 
drawOAvectors(X + 1i*Y, in, p, cm(2,:));

bs_plot = plot(real(bs(:, idx)), imag(bs(:, idx)), 'LineWidth', linewidth);
for i = 1:nlines; cols(i,:) = bs_plot(i).Color; bs_plot(i).Color(4) = linalpha; end
scatter(real(bs(end, idx)), imag(bs(end, idx)), 125, '.k'); 
scatter(real(bs(end,   :)), imag(bs(end,   :)), 125, '.k', 'MarkerEdgeAlpha', linalpha); 
plot(real(ZOA(1:end-10)), imag(ZOA(1:end-10)), 'LineWidth', linewidth, 'color', col);
endline = ZOA(end-1) - ZOA(end);
endpoint = ZOA(end) + 0.06*endline/abs(endline);
plot_arrow(real(endpoint), imag(endpoint), real(ZOA(end)), imag(ZOA(end)),'linewidth', 2, ...
    'color', col,'facecolor', col,'edgecolor', col, 'headwidth',0.5,'headheight',2);

phasespaceplot();
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

% Show the full space:
subplot(3, 3, 9); hold on; box on; axis square;
xlim([-1, 1]); ylim([-1, 1])  
drawOAvectors(X + 1i*Y, in, p, cm(2,:));

bs_plot = plot(real(bs(:, idx)), imag(bs(:, idx)), 'LineWidth', linewidth_total);
for i = 1:nlines; cols(i,:) = bs_plot(i).Color; bs_plot(i).Color(4) = linalpha; end
scatter(real(bs(1, idx)), imag(bs(1, idx)), 100, cols, 'x', 'LineWidth', linewidth_total, 'MarkerEdgeAlpha', linalpha);
bsend_scatter = scatter(real(bs(end, idx)), imag(bs(end, idx)), 100, '.k'); 
scatter(real(bs(end,   :)), imag(bs(end,   :)), 100, '.k', 'MarkerEdgeAlpha', linalpha); 
scatter(real(Zstart), imag(Zstart), 200, [0,0,0], 'x', 'LineWidth', linewidth_total);
scatter(real(ZOA(1)), imag(ZOA(1)), 200, col, '+', 'LineWidth', linewidth_total);
plot(real(ZOA), imag(ZOA), 'LineWidth', linewidth_total, 'color', col);
endline = ZOA(end-2) - ZOA(end);
endpoint = ZOA(end) + 0.08*endline/abs(endline);
plot_arrow(real(endpoint), imag(endpoint), real(ZOA(end)), imag(ZOA(end)),'linewidth', 2, ...
    'color', col,'facecolor', col,'edgecolor', col, 'headwidth',0.7,'headheight',3);

phasespaceplot();
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
disp("Third column")

pause(0.1);
l = legend([Zstart_scatter, ZOAstart_scatter, ZOA_plot, zoa_est, bsstart_scatter(1), bs_plot(1), bsend_scatter], ...
    {'$$Z(0)$$', '$$\bar{Z}(0)$$', '$$\bar{Z}(t)$$', '$$\hat{z}($$\boldmath$$k$$\unboldmath$$, 0)$$',...
    '$$z($$\boldmath$$k$$\unboldmath$$, 0)$$', '$$z($$\boldmath$$k$$\unboldmath$$, t)$$', ...
    '$$z($$\boldmath$$k$$\unboldmath$$, t_{\rm end})$$'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', titlefont);


f1P = [0.5 0.01 0.03 0.03];
set(l,'Position', f1P,'Units', 'normalized');

% Export the figure:
set(findall(gcf,'-property','FontName'),'FontName','Avenir')

exportgraphics(f_mappings,'../Figures/PhaseSpace/Mappings.pdf')

