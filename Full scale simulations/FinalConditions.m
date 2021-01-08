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

%% Theta neurons parameters:
pars.N = 7000;
pars.a_n = 0.666666666666666666667;
tnow = 0; tend = 100;
seed = 3; rng(seed);

%% Draw the problems with mappings, use PSR state:
labelfont = 17;
pars.eta0 = -0.2; pars.delta = 0.1; pars.K = -2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
p = prepareOAparameters(make_scalefreeparameters(pars, 3));

f_finalcond = figure('Renderer', 'painters', 'Position', [50 800 1000 200]); 
Zstart = -0.2 + 1i*0.8; col = p.colorvec;

% With the simple same IC
OAIC = ones(p.Mk,1)*Zstart;
[~, ZOA, bs] = OA_simulatenetwork(0, tend, OAIC, p, odeoptions);

figure; hold on; box on; axis square;
plot(real(bs(end,:)), imag(bs(end,:)))
scatter(real(bs(end,:)), imag(bs(end,:)))

%% Make the histogram on the line:
figure; hold on; box on; axis square;
[edges, counts] = projectToHist(bs(end,:));
xlim([0, 1]);
histogram('BinEdges',edges,'BinCounts',counts, 'Normalization', 'pdf')

% Export the figure:
% set(findall(gcf,'-property','FontName'),'FontName','Avenir')
% exportgraphics(f_finalcond,'../Figures/PhaseSpace/Mappings.pdf')

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
