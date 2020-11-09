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
%     d = gpuDevice(gpuDeviceCount-1);
    d = gpuDevice(1);
    disp(d)
end
initarray = make_GPUhandle();

%% Theta model parameters:
tnow = 0; tend = 50;
h = 0.001;

pars.N = 10000;
pars.a_n = 0.666666666666666666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;

seed = 2; rng(seed);
IC = wrapToPi(randn(pars.N, 1)*1.4);

pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
odeoptions = odeset('RelTol', 1.0e-12,'AbsTol', 1.0e-12);
optimopts = optimoptions('fsolve', 'Display','off', 'Algorithm', 'Levenberg-Marquardt');

% %% 0. Perform a full scale simulation of a FULLY CONNECTED network:
% % The simple DOPRI integration: 
% fdpars = make_fixeddegreeparameters(pars, pars.N);
% [t, thetas] = DOPRI_threshold(@thetaneurons, tnow, tend, IC, h, pars);
% z = orderparameter(thetas);
% disp('Small scale test done')
% 
% % The full scale simulation using the adjacency matrix:
% [tfull, thetasfull] = DOPRI_simulatenetwork(tnow,tend,IC,h,fdpars);
% zfull = orderparameter(thetasfull);
% disp('Full scale test done')
% 
% % The mean field theory for fixed degree networks:
% [T, Z] = ode45(@(t,x) MFR2(t,x,pars), [tnow, tend], gather(zfull(1)), odeoptions);
% disp('Mean field test done')
% 
% % The OA mean field theory:
% fdpars = prepareOAparameters(fdpars);
% z0 = map_thetatozoa(gather(thetasfull(:,1)), fdpars);
% [TOA, ZOA] = OA_simulatenetwork(tnow, tend, gather(z0), fdpars, odeoptions);
% disp('OA mean field test done')
% 
% %% Plotting the results:
% f_fullyconnected = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;
% 
% xlim([tnow, tend]); ylim([0, 1])
% plot(t, abs(z), '-', 'LineWidth', 5, 'Color', '#EDB120');
% plot(tfull, abs(zfull), '-', 'LineWidth', 3, 'Color', '#0072BD');
% plot(T, abs(Z), '-', 'LineWidth', 3, 'Color', '#D95319');
% plot(TOA, abs(ZOA), '-', 'LineWidth', 2, 'Color', '#000000');
% xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
% ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)
% 
% title(sprintf('\\bf Fully connected network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fdpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')
% legend('$$Z(t)_{simple}$$', '$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest', 'Orientation','horizontal')
% exportpdf(f_fullyconnected, '../Figures/InspectMeanFieldFullyConnected.pdf', export);
% close(f_fullyconnected)
% 
% disp('Made fully connected network figure')
% 
% %% 1. Perform a full scale simulation of a fixed degree network:
% netdegree = round(pars.N*0.5);
% 
% % The full scale simulation using the adjacency matrix:
% fdpars = make_fixeddegreeparameters(pars, netdegree);
% [tfull, thetasfull] = DOPRI_simulatenetwork(tnow,tend,IC,h,fdpars);
% zfull = orderparameter(thetasfull);
% disp('Full scale test done')
% 
% % The mean field theory for fixed degree networks:
% [T, Z] = ode45(@(t,x) MFR(t,x,pars), [tnow, tend], gather(zfull(1)), odeoptions);
% disp('Mean field test done')
% 
% % The OA mean field theory:
% fdpars = prepareOAparameters(fdpars);
% z0 = map_thetatozoa(gather(thetasfull(:,1)), fdpars);
% [TOA, ZOA] = OA_simulatenetwork(tnow, tend, z0, fdpars, odeoptions);
% disp('OA mean field test done')
% 
% %% Plotting the results:
% f_fixeddegree = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;
% 
% xlim([tnow, tend]); ylim([0, 1])
% plot(tfull, abs(zfull), '-', 'LineWidth', 4, 'Color', '#0072BD');
% plot(T, abs(Z), '-', 'LineWidth', 3, 'Color', '#D95319');
% plot(TOA, abs(ZOA), '-', 'LineWidth', 2, 'Color', '#000000');
% xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
% ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)
% 
% title(sprintf('\\bf Fixed degree network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fdpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')
% 
% legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest', 'Orientation','horizontal')
% exportpdf(f_fixeddegree, '../Figures/InspectMeanFieldFixedDegree.pdf', export);
% close(f_fixeddegree)
% 
% disp('Made fixed degree network figure')
% 
% %% 2. Perform a full scale simulation of a random network:
% % The full scale simulation using the adjacency matrix:
% netp = 0.3;
% rdpars = make_randomparameters(pars, netp);
% [tfull, thetasfull] = DOPRI_simulatenetwork(tnow,tend,IC,h,rdpars);
% zfull = orderparameter(thetasfull);
% disp('Full scale test done')
% 
% % The OA mean field theory:
% rdpars = prepareOAparameters(rdpars);
% z0 = map_thetatozoa(gather(thetasfull(:,1)), rdpars);
% [TOA, ZOA] = OA_simulatenetwork(tnow, tend, z0, rdpars, odeoptions);
% disp('OA mean field test done')
% 
% %% Plotting the results:
% f_random = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;
% 
% xlim([tnow, tend]); ylim([0, 1])
% plot(tfull, abs(zfull), '-', 'LineWidth', 3, 'Color', '#0072BD');
% plot(TOA, abs(ZOA), '-', 'LineWidth', 2, 'Color', '#000000');
% xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
% ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)
% 
% title(sprintf('\\bf Random network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$p$$ = %0.1f', pars.N, rdpars.meandegree, rdpars.netp), 'FontSize', titlefont, 'Interpreter', 'latex')
% legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest', 'Orientation','horizontal')
% exportpdf(f_random, '../Figures/InspectMeanFieldRandom.pdf', export);
% close(f_random)
% 
% disp('Made random network figure')

%% 3. Perform a full scale simulation of a scale-free network:
degree = 0.3;
IC = wrapToPi(randn(pars.N, 1)*1.2);

% The full scale simulation using the adjacency matrix:
sfpars = make_scalefreeparameters(pars, degree);
[tfull, thetasfull] = DOPRI_simulatenetwork(tnow,tend,IC,h,sfpars);
zfull = orderparameter(thetasfull);
disp('Full scale test done')

% The OA mean field theory:
pars.N = 1000; pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
sfpars = make_scalefreeparameters(pars, degree);
sfpars = prepareOAparameters(sfpars);
% z0 = map_thetatozoa(gather(thetasfull(:,1)), sfpars);
% z0 = orderparameter(IC)*ones(sfpars.Mk,1);
% z0 = ones(sfpars.Mk,1)*(-0.73*1i);
% [TOA, ZOA] = OA_simulatenetwork(tnow, tend, z0, sfpars, odeoptions);
% Limit cycle:
[T, ZOA] = OA_simulatenetwork(0, 100, ones(sfpars.Mk,1)*(-0.2*1i), sfpars);
% plot(real(ZOA), imag(ZOA), 'LineWidth', 2, 'color', col);
ZOA = flip(ZOA(round(numel(T)*0.9):end,:));
[~, pksloc] = findpeaks(abs(ZOA),'MinPeakDistance',100);
idx = pksloc(1):pksloc(3);
plot(real(ZOA(idx)), imag(ZOA(idx)), '-', 'LineWidth', 2);
plot_arrow(real(ZOA(end)), imag(ZOA(end)), real(ZOA(end-2)), imag(ZOA(end-2)),'linewidth', 2, 'headwidth',0.7,'headheight',3);

% disp('OA mean field test done')

close all
f_scalefree = figure('Renderer', 'painters'); hold on; box on; axis square;
drawfixeddegreelimitcycle();
phasespaceplot();
plot(real(zfull), imag(zfull), 'LineWidth', 1, 'Color', '#0072BD');
scatter(real(zfull(1)), imag(zfull(1)), 50, 'x');

plot(real(ZOA), imag(ZOA), '-k','LineWidth', 1);
scatter(real(ZOA(1)), imag(ZOA(1)), 50, '+');

print(f_scalefree, '../Figures/testScaleFree.png', '-dpng', '-r300')
close(f_scalefree)

% %% Plotting the results:
% f_scalefree = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;
% 
% xlim([tnow, tend]); ylim([0, 1])
% plot(tfull, abs(zfull), '-', 'LineWidth', 3, 'Color', '#0072BD');
% plot(TOA, abs(ZOA), '-k', 'LineWidth', 2, 'Color', '#000000');
% xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
% ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)
% 
% title(sprintf('\\bf Scale-free network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$\\gamma$$ = %0.1f', pars.N, sfpars.meandegree, sfpars.degree), 'FontSize', titlefont, 'Interpreter', 'latex')
% legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest', 'Orientation','horizontal')
% exportpdf(f_scalefree, '../Figures/InspectMeanFieldScaleFree.pdf', export);
% close(f_scalefree)
% 
% disp('Made scale-free network figure')

% %% 4. Perform a full scale simulation of a lognorm network:
% % The full scale simulation using the adjacency matrix:
% lnpars = make_lognormparameters(pars, 3, 1, round(pars.N/5));
% [tfull, thetasfull] = DOPRI_simulatenetwork(tnow,tend,IC,h,lnpars);
% zfull = orderparameter(thetasfull);
% disp('Full scale test done')
% 
% % The OA mean field theory:
% lnpars = prepareOAparameters(lnpars);
% odeoptions = odeset('RelTol', 1.0e-9,'AbsTol', 1.0e-9);
% [TOA, ZOA] = OA_simulatenetwork(tnow, tend, IC, lnpars, odeoptions);
% disp('OA mean field test done')
% 
% %% Plotting the results:
% f_lognorm = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;
% 
% xlim([tnow, tend]); ylim([0, 1])
% plot(tfull, abs(zfull), '-', 'LineWidth', 3, 'Color', '#0072BD');
% plot(TOA, abs(ZOA), '-k', 'LineWidth', 2, 'Color', '#000000');
% xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
% ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)
% 
% title(sprintf('\\bf Lognorm network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, lnpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')
% legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest', 'Orientation','horizontal')
% exportpdf(f_lognorm, '../Figures/InspectMeanFieldLogNorm.pdf', export);
% close(f_lognorm)
% 
% disp('Made lognorm network figure')
% 
