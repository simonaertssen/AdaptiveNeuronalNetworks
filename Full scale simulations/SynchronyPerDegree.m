clear all; close all; clc;
% Simulate a full scale fixed degree network and test whether the results
% are correct, with respect to the order parameters.

%% Setup
addpath('../Functions');
addpath('../Mean Field Reductions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 15;
labelfont = 13;
export = true;

%% Make a GPU init handle:
if gpuDeviceCount > 0
    d = gpuDevice(gpuDeviceCount-1);
    disp(d)
end
initarray = make_GPUhandle();

%% Theta model parameters:
tnow = 0; tend = 1;
h = 0.001;

pars.N = 100;
pars.a_n = 0.666666666666666666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;

seed = 2; rng(seed);
IC = - pi/2 * ones(pars.N, 1);

pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
etalower5q = pars.eta0 + pars.delta*tan(pi*(0.05 - 1/2));
etaupper5q = pars.eta0 + pars.delta*tan(pi*(0.95 - 1/2));
lowidx = pars.e < etalower5q;
highidx = pars.e > etaupper5q;

odeoptions = odeset('RelTol', 1.0e-8,'AbsTol', 1.0e-8);
optimopts = optimoptions('fsolve', 'Display','off', 'Algorithm', 'Levenberg-Marquardt');

%% The quantiles of the Cauchy distribution
% figure; hold on;
% x = linspace(-10*pars.delta, 10*pars.delta, 1000);
% plot(x, 1 ./ (pi*pars.delta*(1 + power((x - pars.eta0)./pars.delta, 2))));
% histogram(pars.e, 'Normalization', 'pdf', 'BinEdges', x(1:25:end))
% 
% lower5lim = pars.eta0 + pars.delta*tan(pi*(0.05 - 1/2));
% xline(lower5lim);
% lower5 = pars.e(pars.e < lower5lim);
% histogram(lower5, 'Normalization', 'pdf', 'BinEdges', linspace(x(1), lower5lim, 10))
% 
% upper5lim = pars.eta0 + pars.delta*tan(pi*(0.95 - 1/2));
% xline(upper5lim);
% upper5 = pars.e(pars.e > upper5lim);
% % upper5 = upper5 / sum(upper5);
% histogram(upper5, 'Normalization', 'pdf', 'BinEdges', linspace(upper5lim, x(end), 10))

%% 1. Perform a full scale simulation of a fixed degree network:
netdegree = round(pars.N*0.3);

% The full scale simulation using the adjacency matrix:
fdpars = make_fixeddegreeparameters(pars, netdegree);
[tfull, thetasfull] = DOPRI_simulatenetwork(tnow,tend,IC,h,fdpars);
zfull = orderparameter(thetasfull);
disp('Full scale test done')

% The OA mean field theory:
fdpars = prepareOAparameters(fdpars);
[TOA, ZOA] = OA_simulatenetwork(tnow, tend, IC, fdpars, odeoptions);
disp('OA mean field test done')

% Order parameter per quantile:
zfull_lo = orderparameter(thetasfull(lowidx,:));
zfull_hi = orderparameter(thetasfull(highidx,:));

%% Plotting the results:
f_fixeddegree = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(tfull, abs(zfull), '-', 'LineWidth', 4, 'Color', '#0072BD');
plot(tfull, abs(zfull_lo), '-', 'LineWidth', 2, 'Color', '#D95319');
plot(tfull, abs(zfull_hi), '-', 'LineWidth', 2, 'Color', '#EDB120');
plot(TOA, abs(ZOA), '-', 'LineWidth', 3, 'Color', '#000000');

xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Fixed degree network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fdpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')

legend('$$Z(t)_{A_{ij}}$$', '$$Z(t)_{\eta_{\rm  < 0.05}}$$', '$$Z(t)_{\eta_{\rm > 0.95}}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest', 'Orientation','horizontal')
exportpdf(f_fixeddegree, '../Figures/SynchronyFixedDegree.pdf', export);
close(f_fixeddegree)

disp('Made fixed degree network figure')

%% 2. Perform a full scale simulation of a random network:
% The full scale simulation using the adjacency matrix:
netp = 0.2;
rdpars = make_randomparameters(pars, netp);
[tfull, thetasfull] = DOPRI_simulatenetwork(tnow,tend,IC,h,rdpars);
zfull = orderparameter(thetasfull);
disp('Full scale test done')

% The OA mean field theory:
rdpars = prepareOAparameters(rdpars);
[TOA, ZOA] = OA_simulatenetwork(tnow, tend, IC, rdpars, odeoptions);
disp('OA mean field test done')

% Order parameter per quantile:
zfull_lo = orderparameter(thetasfull(lowidx));
zfull_hi = orderparameter(thetasfull(highidx));

% Order parameter per degree:
nbins = 4;
[N,edges] = histcounts(rdpars.degrees_in, nbins);
zdegrees = zeros(nbins, numel(tfull));  
for i = 1:nbins
    idx = rdpars.degrees_in > edges(i) & rdpars.degrees_in < edges(i+1);
    zdegrees(i,:) = orderparameter(thetasfull(idx,:));
end

%% Plotting the results:
f_random = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(tfull, abs(zfull), '-', 'LineWidth', 3, 'Color', '#0072BD');
plot(tfull, abs(zfull_lo), '-', 'LineWidth', 2, 'Color', '#D95319');
plot(tfull, abs(zfull_hi), '-', 'LineWidth', 2, 'Color', '#EDB120');
plot(TOA, abs(ZOA), '-', 'LineWidth', 2, 'Color', '#000000');

for i = 1:nbins
    plot(tfull, abs(zdegrees(i,:)), '-', 'LineWidth', 2, 'Color', [0, 1/nbins*i, 0])
end

xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Random network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$p$$ = %0.1f', pars.N, rdpars.meandegree, rdpars.netp), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{A_{ij}}$$', '$$Z(t)_{\eta_{\rm  < 0.05}}$$', '$$Z(t)_{\eta_{\rm > 0.95}}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest', 'Orientation','horizontal')
exportpdf(f_random, '../Figures/SynchronyRandom.pdf', export);
% close(f_random)

disp('Made random network figure')

%% 3. Perform a full scale simulation of a scale-free network:
% The full scale simulation using the adjacency matrix:
degree = 3;

sfpars = make_scalefreeparameters(pars, degree);
[tfull, thetasfull] = DOPRI_simulatenetwork(tnow,tend,IC,h,sfpars);
zfull = orderparameter(thetasfull);
disp('Full scale test done')

% The OA mean field theory:
sfpars = prepareOAparameters(sfpars);
[TOA, ZOA] = OA_simulatenetwork(tnow, tend, IC, sfpars, odeoptions);
disp('OA mean field test done')

% Order parameter per quantile:
zfull_lo = orderparameter(thetasfull(lowidx));
zfull_hi = orderparameter(thetasfull(highidx));

% Order parameter per degree:
nbins = 4;
[N,edges] = histcounts(rdpars.degrees_in, nbins);
zdegrees = zeros(nbins, numel(tfull));  
for i = 1:nbins
    idx = rdpars.degrees_in > edges(i) & rdpars.degrees_in < edges(i+1);
    zdegrees(i,:) = orderparameter(thetasfull(idx,:));
end


%% Plotting the results:
f_scalefree = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(tfull, abs(zfull), '-', 'LineWidth', 3, 'Color', '#0072BD');
plot(tfull, abs(zfull_lo), '-', 'LineWidth', 2, 'Color', '#D95319');
plot(tfull, abs(zfull_hi), '-', 'LineWidth', 2, 'Color', '#EDB120');
plot(TOA, abs(ZOA), '-', 'LineWidth', 2, 'Color', '#000000');

for i = 1:nbins
    plot(tfull, abs(zdegrees(i,:)), '-', 'LineWidth', 2, 'Color', [0, 1/nbins*i, 0])
end

title(sprintf('\\bf Scale-free network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$\\gamma$$ = %0.1f', pars.N, sfpars.meandegree, sfpars.degree), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'northwest', 'Orientation','horizontal')
exportpdf(f_scalefree, '../Figures/SynchronyScaleFree.pdf', export);
close(f_scalefree)

disp('Made scale-free network figure')

