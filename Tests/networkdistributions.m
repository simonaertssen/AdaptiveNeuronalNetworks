close all; clear all; clc;
addpath('../Functions');
labelfont = 15;
titlefont = 15;

%% Test the different networks and their statistical properties
% We need to evaluate whether some of the random numbers we are pulling are
% actually following the right pdf.

pars.N = 10000;
meandegreetoget = 100;
f_1Dpdfs = figure('Renderer', 'painters', 'Position', [50 800 900 400]); box on; hold on;

%% A 1D random network: CORRECT
randompars = make_randomparameters(pars, meandegreetoget/(pars.N - 1));

if round(sum(randompars.P(1:pars.N))) == pars.N; disp('Sum is correct'); end
if abs(randompars.meandegree - sum(randompars.degrees_i)/pars.N) < 1.0e-3*randompars.meandegree; disp('Mean degree is correct'); end

limits = round(randompars.meandegree + sqrt(randompars.meandegree)*[-3, 3]);

subplot(2,3,2); hold on; grid on; axis square; box on;
title("Random", 'FontSize', titlefont)
xlim(limits); x = limits(1):limits(2);
histogram(randompars.degrees_i, 'Normalization', 'pdf');
plot(x, randompars.P(x)/pars.N, 'LineWidth', 3);
xline(randompars.meandegree, 'k--', 'LineWidth', 2)
xlabel('Degree', 'Interpreter', 'none', 'FontSize', labelfont)

subplot(2,3,5); hold on; grid on; axis square; box on;
scatter(randompars.degrees_i, randompars.degrees_o, '.k');
xlabel('$$k^{\rm in}$$', 'Interpreter', 'latex', 'FontSize', labelfont)

%% A 1D scale free network: CORRECT
kmin = 68; kmax = 200;
scalefreepars = make_scalefreeparameters(pars, 3, kmin, kmax);
scalefreepars.meandegree

if round(sum(scalefreepars.P(1:pars.N))) == pars.N; disp('Sum is correct'); end
if abs(scalefreepars.meandegree - sum(scalefreepars.degrees_i)/pars.N) < 1.0e-3*scalefreepars.meandegree; disp('Mean degree is correct'); end

limits = [kmin, kmax];

plt2 = subplot(2,3,3); hold on; grid on; axis square; box on;
title("Scale-free", 'FontSize', titlefont)

xlim([limits(1) - 10, limits(2) + 10]); x = limits(1):limits(2);
histogram(scalefreepars.degrees_i, 'Normalization', 'pdf', 'BinEdges', linspace(kmin, kmax, 20));
plot(x, scalefreepars.P(x)/pars.N, 'LineWidth', 3);
xline(scalefreepars.meandegree, 'k--', 'LineWidth', 2)
xlabel('Degree', 'Interpreter', 'none', 'FontSize', labelfont)
pos = plt2.Position; plt2.Position = [pos(1) - 0.09, pos(2), pos(3), pos(4)];

plt2 = subplot(2,3,6); hold on; grid on; axis square; box on;
xlim(limits); ylim(limits);
scatter(scalefreepars.degrees_i, scalefreepars.degrees_o, '.k');
xlabel('$$k^{\rm in}$$', 'Interpreter', 'latex', 'FontSize', labelfont)
pos = plt2.Position; plt2.Position = [pos(1) - 0.09, pos(2), pos(3), pos(4)];

%% A 1D fixed degree network / diracnet: CORRECT
fixeddegreepars = make_fixeddegreeparameters(pars, meandegreetoget);

if sum(fixeddegreepars.P(1:pars.N)) == pars.N; disp('Sum is correct'); end
if fixeddegreepars.meandegree == sum(fixeddegreepars.degrees_i)/pars.N; disp('Mean degree is correct'); end

limits = round(fixeddegreepars.meandegree + fixeddegreepars.meandegree*[-0.02, 0.02]);

plt = subplot(2,3,1); hold on; grid on; axis square; box on;
title("Fixed degree", 'FontSize', titlefont)
xlim(limits); x = limits(1):0.01:limits(2);
histogram(fixeddegreepars.degrees_i, 'Normalization', 'pdf');
plot(x, fixeddegreepars.P(x)/pars.N, 'LineWidth', 3);
xline(fixeddegreepars.meandegree, 'k--', 'LineWidth', 2)
xlabel('Degree', 'Interpreter', 'none', 'FontSize', labelfont)
ylabel('Density [%]', 'Interpreter', 'none', 'FontSize', labelfont);
pos = plt.Position; plt.Position = [pos(1) + 0.09, pos(2), pos(3), pos(4)];

plt = subplot(2,3,4); hold on; grid on; axis square; box on;
xlim([fixeddegreepars.meandegree-1, fixeddegreepars.meandegree+1]);
ylim([fixeddegreepars.meandegree-1, fixeddegreepars.meandegree+1]);
xticks([fixeddegreepars.meandegree-1, fixeddegreepars.meandegree, fixeddegreepars.meandegree+1])
yticks([fixeddegreepars.meandegree-1, fixeddegreepars.meandegree, fixeddegreepars.meandegree+1])
scatter(fixeddegreepars.degrees_i, fixeddegreepars.degrees_o, '.k');
xlabel('$$k^{\rm in}$$', 'Interpreter', 'latex', 'FontSize', labelfont)
ylabel('$$k^{\rm out}$$', 'Interpreter', 'latex', 'FontSize', labelfont)
pos = plt.Position; plt.Position = [pos(1) + 0.09, pos(2), pos(3), pos(4)];

%% Produce 1D figure
print(f_1Dpdfs, '../Figures/Distributions/1D.png', '-dpng', '-r300')

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

%% Now the 2D pdfs:
%% A 2D fixed degree network / diracnet: CORRECT
fixeddegreepars = make_fixeddegreeparameters(pars, meandegreetoget);

fixeddegreepars.P2D = @(x,y) fixeddegreepars.P(x) .* fixeddegreepars.P(y) / pars.N;

if sum(fixeddegreepars.P2D(1:pars.N, 1:pars.N), 'all') == pars.N; disp('Sum is correct'); end

%%
[x,y] = meshgrid(minmax, minmax);
idx = round(linspace(1, numel(minmax), 25));
sum(fixeddegreepars.P2D(1:pars.N, 1:pars.N), 'all')

%%
minmax = linspace(meandegreetoget - 10, meandegreetoget + 10, 21);

histogram2(fixeddegreepars.degrees_i, fixeddegreepars.degrees_o, minmax, minmax, 'Normalization', 'pdf'); % Normal degree vectors from before
