clear all; close all; clc;
% Simulate a full scale fixed degree network and test whether the results
% are correct, with respect to the order parameters.

%% Setup
addpath('../Functions');
addpath('../Mean Field Reductions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 16;
labelfont = 15;
export = true;

%% Make a GPU init handle:
if gpuDeviceCount > 0
    d = gpuDevice(gpuDeviceCount-1);
end
initarray = make_GPUhandle();

%% Theta model parameters:
tnow = 0; tend = 10;
h = 0.001;

p.N = 15000;
p.a_n = 0.666666666666666666667;

seed = 2; rng(seed);

PSRp = p; PSRp.text = '\it PSR';
PSRp.IC = pi*ones(p.N, 1) - pi; 
PSRp.eta0 = -0.9; PSRp.delta = 0.8; PSRp.K = -2;
PSRp.e = randcauchy(seed, PSRp.eta0, PSRp.delta, PSRp.N);

PSSp = p; PSSp.text = '\it PSS';
PSSp.IC = wrapToPi(randn(p.N, 1)*1.9); %pi*ones(p.N, 1) - pi; % wrapToPi(2*pi*rand(p.N, 1)-pi);
PSSp.eta0 = 0.5; PSSp.delta = 0.5; PSSp.K = 1;
PSSp.e = randcauchy(seed, PSSp.eta0, PSSp.delta, PSSp.N);

CPWp = p; CPWp.text = '\it CPW';
CPWp.IC = wrapToPi(randn(p.N, 1)*1.55);
CPWp.eta0 = 10.75; CPWp.delta = 0.5; CPWp.K = -9;
CPWp.e = randcauchy(seed, CPWp.eta0, CPWp.delta, CPWp.N);

odeoptions = odeset('RelTol', 1.0e-9,'AbsTol', 1.0e-9);

%% 0. Perform a full scale simulation of a FULLY CONNECTED network:
% f_fullyconnected = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;
% xlim([tnow, tend]); ylim([0, 1]);
% A = 0;
% 
% for i = 1:3
%     if i == 1
%         pars = PSRp;
%     elseif i == 2
%         pars = PSSp;
%     elseif i == 3
%         pars = CPWp;
%     else
%         warning('no pars?')
%     end
%     
%     % The simple DOPRI integration:
%     fdpars = make_fixeddegreeparameters(pars, pars.N);
% 
%     % The full scale simulation using the adjacency matrix:
%     [tfull, thetasfull, A] = DOPRI_simulatenetwork(tnow,tend,pars.IC,h,fdpars,A);
%     zfull = orderparameter(thetasfull);
%     disp('Full scale test done')
% 
%     % The mean field theory for fixed degree networks:
%     [T, Z] = ode45(@(t,x) MFR2(t,x,pars), [tnow, tend], gather(zfull(1)), odeoptions);
%     disp('Mean field test done')
% 
%     % The OA mean field theory:
%     fdpars = prepareOAparameters(fdpars);
%     z0 = orderparameter(pars.IC)*ones(fdpars.Mk,1);
%     [TOA, ZOA] = OA_simulatenetwork(tnow, tend, gather(z0), fdpars, odeoptions);
%     disp('OA mean field test done')
%     
%     % Plotting the results:
%     plot(tfull, abs(zfull), '-', 'LineWidth', 3, 'Color', '#0072BD');
%     plot(T, abs(Z), '-', 'LineWidth', 2, 'Color', '#D95319');
%     plot(TOA, abs(ZOA), '-', 'LineWidth', 1, 'Color', '#000000');
% end
% 
% text(tend*0.99, 0.99, PSRp.text, 'FontSize', labelfont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
% text(tend*0.99, 0.20, PSSp.text, 'FontSize', labelfont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')
% text(tend*0.99, 0.70, CPWp.text, 'FontSize', labelfont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')
% 
% xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
% ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)
% title(sprintf('\\bf Fully connected network: $$N$$ = %d, $$\\langle k \\rangle$$ = %d', pars.N, fdpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')
% 
% legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southeast', 'Orientation','horizontal')
% exportpdf(f_fullyconnected, '../Figures/InspectMeanFieldFullyConnected.pdf', export);
% close(f_fullyconnected)

% disp('Made fully connected network figure')

%% 1. Perform a full scale simulation of a fixed degree network:
netdegree = round(p.N*0.5);

f_fixeddegree = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;
xlim([tnow, tend]); ylim([0, 1]);
A = 0;

for i = 1:3
    if i == 1
        pars = PSRp;
    elseif i == 2
        pars = PSSp;
    elseif i == 3
        pars = CPWp;
    else
        warning('no pars?')
    end

    % The full scale simulation using the adjacency matrix:
    fdpars = make_fixeddegreeparameters(pars, netdegree);
    [tfull, thetasfull, A] = DOPRI_simulatenetwork(tnow,tend,pars.IC,h,fdpars,A);
    zfull = orderparameter(thetasfull);
    disp('Full scale test done')

    % The mean field theory for fixed degree networks:
    [T, Z] = ode45(@(t,x) MFR(t,x,pars), [tnow, tend], gather(zfull(1)), odeoptions);
    disp('Mean field test done')

    % The OA mean field theory:
    fdpars = prepareOAparameters(fdpars);
    z0 = orderparameter(pars.IC)*ones(fdpars.Mk,1);
    [TOA, ZOA] = OA_simulatenetwork(tnow, tend, z0, fdpars, odeoptions);
    disp('OA mean field test done')

    % Plotting the results:
    plot(tfull, abs(zfull), '-', 'LineWidth', 3, 'Color', '#0072BD');
    plot(T, abs(Z), '-', 'LineWidth', 2, 'Color', '#D95319');
    plot(TOA, abs(ZOA), '-', 'LineWidth', 1, 'Color', '#000000');
end

text(tend*0.99, 0.99, PSRp.text, 'FontSize', labelfont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
text(tend*0.99, 0.20, PSSp.text, 'FontSize', labelfont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')
text(tend*0.99, 0.70, CPWp.text, 'FontSize', labelfont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')

xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Fixed degree network: $$N$$ = %d, $$\\langle k \\rangle$$ = %d', pars.N, fdpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')

legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z}(t)_{MF}$$', '$$\overline{Z}(t)_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southeast', 'Orientation','horizontal')
exportpdf(f_fixeddegree, '../Figures/InspectMeanFieldFixedDegree.pdf', export);
close(f_fixeddegree)

disp('Made fixed degree network figure')

%% 2. Perform a full scale simulation of a random network:
netp = 0.3;

f_random = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;
xlim([tnow, tend]); ylim([0, 1]);
A = 0;

for i = 1:3
    if i == 1
        pars = PSRp;
    elseif i == 2
        pars = PSSp;
    elseif i == 3
        pars = CPWp;
    else
        warning('no pars?')
    end

    % The full scale simulation using the adjacency matrix:
    rdpars = make_randomparameters(pars, netp);
    [tfull, thetasfull, A] = DOPRI_simulatenetwork(tnow,tend,pars.IC,h,rdpars,A);
    zfull = orderparameter(thetasfull);
    disp('Full scale test done')

    % The OA mean field theory:
    rdpars = prepareOAparameters(rdpars);
    z0 = orderparameter(pars.IC)*ones(rdpars.Mk,1);
    [TOA, ZOA] = OA_simulatenetwork(tnow, tend, z0, rdpars, odeoptions);
    disp('OA mean field test done')

    % Plotting the results:
    plot(tfull, abs(zfull), '-', 'LineWidth', 3, 'Color', '#0072BD');
    plot(TOA, abs(ZOA), '-', 'LineWidth', 2, 'Color', '#000000');
end

text(tend*0.99, 0.99, PSRp.text, 'FontSize', labelfont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
text(tend*0.99, 0.20, PSSp.text, 'FontSize', labelfont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')
text(tend*0.99, 0.70, CPWp.text, 'FontSize', labelfont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')

xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Random network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$p$$ = %0.1f', pars.N, rdpars.meandegree, rdpars.netp), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z}(t)_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southeast', 'Orientation','horizontal')
exportpdf(f_random, '../Figures/InspectMeanFieldRandom.pdf', export);
close(f_random)

disp('Made random network figure')

%% 3. Perform a full scale simulation of a scale-free network:
degree = 2.04;
% IC = wrapToPi(randn(pars.N, 1)*1.2);

f_scalefree = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;
xlim([tnow, tend]); ylim([0, 1]);
A = 0;

for i = 1:3
    if i == 1
        pars = PSRp;
    elseif i == 2
        pars = PSSp;
    elseif i == 3
        pars = CPWp;
    else
        warning('no pars?')
    end

    % The full scale simulation using the adjacency matrix:
    sfpars = make_scalefreeparameters(pars, degree);
    [tfull, thetasfull, A] = DOPRI_simulatenetwork(tnow,tend,pars.IC,h,sfpars,A);
    zfull = orderparameter(thetasfull);
    disp('Full scale test done')

    % The OA mean field theory:
    sfpars = prepareOAparameters(sfpars);
    z0 = orderparameter(pars.IC)*ones(sfpars.Mk,1);
    [TOA, ZOA] = OA_simulatenetwork(tnow, tend, z0, sfpars, odeoptions);
    disp('OA mean field test done')

    % Plotting the results:
    plot(tfull, abs(zfull), '-', 'LineWidth', 3, 'Color', '#0072BD');
    plot(TOA, abs(ZOA), '-k', 'LineWidth', 2, 'Color', '#000000');
end

text(tend*0.99, 0.99, PSRp.text, 'FontSize', labelfont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
text(tend*0.99, 0.20, PSSp.text, 'FontSize', labelfont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')
text(tend*0.99, 0.70, CPWp.text, 'FontSize', labelfont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')

xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Scale-free network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$\\gamma$$ = %0.2f', pars.N, sfpars.meandegree, sfpars.degree), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z}(t)_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southeast', 'Orientation','horizontal')
exportpdf(f_scalefree, '../Figures/InspectMeanFieldScaleFree.pdf', export);
close(f_scalefree)

disp('Made scale-free network figure')

%% 4. Perform a full scale simulation of a lognorm network:
f_lognorm = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;
xlim([tnow, tend]); ylim([0, 1]);
A = 0;

for i = 1:3
    if i == 1
        pars = PSRp;
    elseif i == 2
        pars = PSSp;
    elseif i == 3
        pars = CPWp;
    else
        warning('no pars?')
    end

    % The full scale simulation using the adjacency matrix:
    mu = 3; sigma = 1; kmin = round(pars.N/5);
    lnpars = make_lognormparameters(pars, mu, sigma, kmin);
    [tfull, thetasfull, A] = DOPRI_simulatenetwork(tnow,tend,pars.IC,h,lnpars,A);
    zfull = orderparameter(thetasfull);
    disp('Full scale test done')

    % The OA mean field theory:
    lnpars = prepareOAparameters(lnpars);
    z0 = orderparameter(pars.IC)*ones(lnpars.Mk,1);
    [TOA, ZOA] = OA_simulatenetwork(tnow, tend, z0, lnpars, odeoptions);
    disp('OA mean field test done')

    % Plotting the results:
    plot(tfull, abs(zfull), '-', 'LineWidth', 3, 'Color', '#0072BD');
    plot(TOA, abs(ZOA), '-k', 'LineWidth', 2, 'Color', '#000000');
end

text(tend*0.99, 0.99, PSRp.text, 'FontSize', labelfont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
text(tend*0.99, 0.20, PSSp.text, 'FontSize', labelfont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')
text(tend*0.99, 0.70, CPWp.text, 'FontSize', labelfont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')

xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Lognorm network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$\\mu$$ = %0.1f, $$\\sigma$$ = %0.1f, $$k_{\\rm min}$$ = %d', pars.N, lnpars.meandegree, mu, sigma, kmin), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z}(t)_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southeast', 'Orientation','horizontal')
exportpdf(f_lognorm, '../Figures/InspectMeanFieldLogNorm.pdf', export);
close(f_lognorm)

disp('Made lognorm network figure')

