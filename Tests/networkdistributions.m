close all; clear all; clc;
addpath('../Functions');
labelfont = 15;
titlefont = 15;

make1Dplots = true;
make2Dplots = true;

%% Test the different networks and their statistical properties
% We need to evaluate whether some of the random numbers we are pulling are
% actually following the right pdf.

pars.N = 10000;

meandegreetoget = 100;

%% Figure handle
if make1Dplots == true
f_1Dpdfs = figure('Renderer', 'painters', 'Position', [50 800 1000 200]); box on; hold on;

%% A 1D random network: CORRECT
randompars = make_randomparameters(pars, meandegreetoget/(pars.N - 1));

if round(sum(randompars.P(1:pars.N))) == pars.N; disp('Sum is correct'); end
if abs(randompars.meandegree - sum(randompars.degrees_i)/pars.N) < 1.0e-3*randompars.meandegree; disp('Mean degree is correct'); end

limits = round(randompars.meandegree + sqrt(randompars.meandegree)*[-3, 3]);

subplot(1,3,2); hold on; grid on; axis square; box on;
title("Random", 'FontSize', titlefont)
xlim(limits); x = limits(1):limits(2);
histogram(randompars.degrees_i, 'Normalization', 'pdf');
plot(x, randompars.P(x)/pars.N, 'LineWidth', 3);
xline(randompars.meandegree, 'k--', 'LineWidth', 2)
xlabel('Degree', 'Interpreter', 'none', 'FontSize', labelfont)

% subplot(1,3,5); hold on; grid on; axis square; box on;
% scatter(randompars.degrees_i, randompars.degrees_o, '.k');
% xlabel('$$k^{\rm in}$$', 'Interpreter', 'latex', 'FontSize', labelfont)

%% A 1D scale free network: CORRECT
kmin = 68; kmax = 200;
scalefreepars = make_scalefreeparameters(pars, 3, kmin, kmax);
scalefreepars.meandegree

if round(sum(scalefreepars.P(1:pars.N))) == pars.N; disp('Sum is correct'); end
if abs(scalefreepars.meandegree - sum(scalefreepars.degrees_i)/pars.N) < 1.0e-3*scalefreepars.meandegree; disp('Mean degree is correct'); end

limits = [kmin, kmax];

plt2 = subplot(1,3,3); hold on; grid on; axis square; box on;
title("Scale-free", 'FontSize', titlefont)

xlim([limits(1) - 10, limits(2) + 10]); x = limits(1):limits(2);
histogram(scalefreepars.degrees_i, 'Normalization', 'pdf', 'BinEdges', linspace(kmin, kmax, 20));
plot(x, scalefreepars.P(x)/pars.N, 'LineWidth', 3);
xline(scalefreepars.meandegree, 'k--', 'LineWidth', 2)
xlabel('Degree', 'Interpreter', 'none', 'FontSize', labelfont)
pos = plt2.Position; plt2.Position = [pos(1) - 0.09, pos(2), pos(3), pos(4)];

% plt2 = subplot(2,3,6); hold on; grid on; axis square; box on;
% xlim(limits); ylim(limits);
% scatter(scalefreepars.degrees_i, scalefreepars.degrees_o, '.k');
% xlabel('$$k^{\rm in}$$', 'Interpreter', 'latex', 'FontSize', labelfont)
% pos = plt2.Position; plt2.Position = [pos(1) - 0.09, pos(2), pos(3), pos(4)];

%% A 1D fixed degree network / diracnet: CORRECT
fixeddegreepars = make_fixeddegreeparameters(pars, meandegreetoget);

if sum(fixeddegreepars.P(1:pars.N)) == pars.N; disp('Sum is correct'); end
if fixeddegreepars.meandegree == sum(fixeddegreepars.degrees_i)/pars.N; disp('Mean degree is correct'); end

limits = round(fixeddegreepars.meandegree + fixeddegreepars.meandegree*[-0.02, 0.02]);

plt = subplot(1,3,1); hold on; grid on; axis square; box on;
title("Fixed degree", 'FontSize', titlefont)
xlim(limits); x = limits(1):0.01:limits(2);
histogram(fixeddegreepars.degrees_i, 'Normalization', 'pdf');
plot(x, fixeddegreepars.P(x)/pars.N, 'LineWidth', 3);
xline(fixeddegreepars.meandegree, 'k--', 'LineWidth', 2)
xlabel('Degree', 'Interpreter', 'none', 'FontSize', labelfont)
ylabel('Density [%]', 'Interpreter', 'none', 'FontSize', labelfont);
pos = plt.Position; plt.Position = [pos(1) + 0.09, pos(2), pos(3), pos(4)];

% plt = subplot(2,3,4); hold on; grid on; axis square; box on;
% xlim([fixeddegreepars.meandegree-1, fixeddegreepars.meandegree+1]);
% ylim([fixeddegreepars.meandegree-1, fixeddegreepars.meandegree+1]);
% xticks([fixeddegreepars.meandegree-1, fixeddegreepars.meandegree, fixeddegreepars.meandegree+1])
% yticks([fixeddegreepars.meandegree-1, fixeddegreepars.meandegree, fixeddegreepars.meandegree+1])
% scatter(fixeddegreepars.degrees_i, fixeddegreepars.degrees_o, '.k');
% xlabel('$$k^{\rm in}$$', 'Interpreter', 'latex', 'FontSize', labelfont)
% ylabel('$$k^{\rm out}$$', 'Interpreter', 'latex', 'FontSize', labelfont)
% pos = plt.Position; plt.Position = [pos(1) + 0.09, pos(2), pos(3), pos(4)];

%% Test: a 1D lognormal network - inspired by scale-free CORRECT
% lognormpars = make_lognormparameters(pars, 3, 1, 500);
% 
% if round(sum(lognormpars.P(1:pars.N))) == pars.N; disp('Sum is correct'); end
% if abs(lognormpars.meandegree - sum(lognormpars.degrees_i)/pars.N) < 1.0e-3*lognormpars.meandegree; disp('Mean degree is correct'); end
% 
% limits = [lognormpars.kmin-10, lognormpars.kmin+200];
% 
% figure; hold on; grid on;
% x = limits(1):0.01:limits(2);
% xlim([limits(1), limits(2)]);
% histogram(lognormpars.degrees_i, 'Normalization', 'pdf');
% plot(x, lognormpars.P(x)/pars.N);
% xline(lognormpars.meandegree, 'k--')

%% Produce 1D figure
set(findall(gcf,'-property','FontName'),'FontName','Avenir')
exportgraphics(f_1Dpdfs,'../Figures/Distributions/1D.pdf')

close(f_1Dpdfs)
end
%% Now the 2D pdfs:
meandegreetoget = pars.N/5;

if make2Dplots == true
f_2Dpdfs = figure('Renderer', 'painters', 'Position', [50 800 900 400]); box on; hold on;

%% A 2D fixed degree network / diracnet
fixeddegreepars = make_fixeddegreeparameters(pars, meandegreetoget);

fixeddegreepars.P2D = @(x,y) fixeddegreepars.P(x) .* fixeddegreepars.P(y) / pars.N;

if sum(fixeddegreepars.P2D(1:pars.N, 1:pars.N), 'all') == pars.N; disp('Sum is correct'); end

%%
subplot(1,3,1); hold on; grid on; box on;
title('Fixed-degree', 'FontSize', titlefont);
xlabel('\boldmath$k^{\rm in}$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('\boldmath$k^{\rm out}$', 'Interpreter', 'latex', 'FontSize', labelfont);
zlabel('Density', 'FontSize', labelfont, 'HorizontalAlignment', 'left');
view(110, 10)

lo = meandegreetoget - 5;
hi = meandegreetoget + 5;

% The 3D histogram
minmax = linspace(lo, hi, hi - lo);
xlim([minmax(1), minmax(end)])
ylim([minmax(1), minmax(end)])
hist = histogram2(fixeddegreepars.degrees_i, fixeddegreepars.degrees_o, minmax, minmax, 'Normalization', 'pdf','ShowEmptyBins','on'); % Normal degree vectors from before
hist.FaceAlpha = 0.75;
colormap(jet)

% The pdf
plot3([meandegreetoget, meandegreetoget], [meandegreetoget, meandegreetoget], [0, max(hist.Values, [], 'all')], 'k', 'LineWidth', 2)
colormap(jet)

% The 2D scatter plot underneath
binwidth = hist.XBinEdges(2) - hist.XBinEdges(1);
binsvalues = zeros(numel(minmax), numel(minmax));
binsvalues(1:numel(minmax)-1, 1:numel(minmax)-1) = hist.Values;
h = pcolor(minmax,minmax,binsvalues);
h.ZData = -0.5*max(binsvalues, [], 'all')*ones(size(binsvalues));

ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
colormap(jet)

%% A 2D random network:
randompars = make_randomparameters(pars, meandegreetoget/(pars.N - 1));

% randompars.P2D = @(x,y) randompars.P(x) .* randompars.P(y) / sum(randompars.P(x).*randompars.P(y), 'all') * pars.N;
% Pnorm = randompars.P2D(1:pars.N,1:pars.N);
% randompars.P2D = @(x,y) randompars.P2D(x,y)/sum(Pnorm, 'all') * pars.N;

randompars.P2D = @(x,y) randompars.P(x).*randompars.P(y);
Pnorm = sum(randompars.P2D(1:pars.N,1:pars.N), 'all')/pars.N; 
randompars.P2D = @(x,y) randompars.P2D(x,y)/Pnorm;

sum(randompars.P2D(1:pars.N, 1:pars.N), 'all')
if abs(sum(randompars.P2D(1:pars.N, 1:pars.N), 'all') - pars.N) < 1.e-3; disp('Sum is correct'); end

%%
subplot(1,3,2); hold on; grid on; box on;
% figure; hold on;

title('Random', 'FontSize', titlefont);
xlabel('\boldmath$k^{\rm in}$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('\boldmath$k^{\rm out}$', 'Interpreter', 'latex', 'FontSize', labelfont);
% zlabel('Density', 'FontSize', labelfont, 'HorizontalAlignment', 'left');
view(110, 10)

% The 3D histogram
minmax = linspace(randompars.meandegree - round(4*sqrt(randompars.meandegree)), randompars.meandegree + round(4*sqrt(randompars.meandegree)), 16);
xlim([minmax(1), minmax(end)])
ylim([minmax(1), minmax(end)])
hist = histogram2(randompars.degrees_i, randompars.degrees_o, minmax, minmax, 'Normalization', 'pdf','FaceAlpha',1,'ShowEmptyBins','on'); % Normal degree vectors from before
colormap(jet)

% The pdf
vec = linspace(minmax(1), minmax(end), minmax(end) - minmax(1) + 1);
[x,y] = meshgrid(vec, vec);
idx = round(linspace(1, numel(vec), 25));
surf(x(idx,idx),y(idx,idx),randompars.P2D(x(idx,idx),y(idx,idx))/pars.N^2*Pnorm,'FaceAlpha',0.5,'EdgeColor','none');

colormap(jet)

% The 2D scatter plot underneath
binwidth = hist.XBinEdges(2) - hist.XBinEdges(1);
binsvalues = zeros(numel(minmax), numel(minmax));
binsvalues(1:numel(minmax)-1, 1:numel(minmax)-1) = hist.Values;
h = pcolor(minmax,minmax,binsvalues);
h.ZData = -0.5*max(binsvalues, [], 'all')*ones(size(binsvalues));

ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
colormap(jet)

%% A 2D scalefree network 
scalefreepars = make_scalefreeparameters(pars, 4.305);
scalefreepars.meandegree

% scalefreepars.P2D = @(x,y) scalefreepars.P(x).*scalefreepars.P(y);
% Pnorm = sum(scalefreepars.P2D(1:pars.N,1:pars.N), 'all')/pars.N; 
% scalefreepars.P2D = @(x,y) scalefreepars.P2D(x,y)/Pnorm;

scalefreepars.P2D = @(x,y) P2D(x, y, scalefreepars.kmin, scalefreepars.kmax, scalefreepars.degree, scalefreepars.N);

vec = linspace(scalefreepars.kmin, scalefreepars.kmax, scalefreepars.kmax-scalefreepars.kmin+1);
[x,y] = meshgrid(vec, vec);
sum(scalefreepars.P2D(x,y), 'all') 
if abs(sum(scalefreepars.P2D(x,y), 'all') - pars.N) < 1.e-3; disp('Sum is correct'); end

%%
% subplot(1,3,3); hold on; grid on; box on;
figure; hold on; box on;
title('Scale-free', 'FontSize', titlefont);
xlabel('\boldmath$k^{\rm in}$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('\boldmath$k^{\rm out}$', 'Interpreter', 'latex', 'FontSize', labelfont);
% zlabel('Density', 'FontSize', labelfont, 'HorizontalAlignment', 'left');
view(110, 10)

% The 3D histogram
minmax = linspace(scalefreepars.kmin, scalefreepars.kmax, 16);
hist = histogram2(scalefreepars.degrees_i, scalefreepars.degrees_o, minmax, minmax, 'Normalization', 'pdf'); % Normal degree vectors from before
xlim([minmax(1), minmax(end)])
ylim([minmax(1), minmax(end)])
colormap(jet)

% The pdf
minmaxpdf = linspace(scalefreepars.kmin, scalefreepars.kmax, scalefreepars.kmax-scalefreepars.kmin+1);
[x,y] = meshgrid(minmax, minmax);
idx = round(linspace(1, numel(minmax), 25));
surf(x(idx,idx),y(idx,idx),scalefreepars.P2D(x(idx,idx),y(idx,idx))/pars.N,'FaceAlpha',0.5,'EdgeColor','none');
colormap(jet)

% The 2D scatter plot underneath
binwidth = hist.XBinEdges(2) - hist.XBinEdges(1);
binsvalues = zeros(numel(minmax), numel(minmax));
binsvalues(1:numel(minmax)-1, 1:numel(minmax)-1) = hist.Values;
h = pcolor(minmax,minmax,binsvalues);
h.ZData = -0.5*max(binsvalues, [], 'all')*ones(size(binsvalues));

ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
colormap(jet)

%% Produce 2D figure
set(findall(gcf,'-property','FontName'),'FontName','Avenir')
exportgraphics(f_2Dpdfs,'../Figures/Distributions/2D.pdf')

close(f_2Dpdfs)
end
function P = P2D(X,Y, kmin, kmax, degree, N)
    P = X.^(-degree) .* Y.^(-degree);
    
    vec = linspace(kmin, kmax, kmax-kmin+1);
    [x,y] = meshgrid(vec, vec);
    Pnorm = x.^(-degree) .* y.^(-degree);
    P = P/sum(Pnorm, 'all')*N;
end

