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
tnow = 0; tend = 20;
h = 0.005;

pars.N = 10000;
pars.a_n = 0.666666666666666666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;

seed = 2; rng(seed);
IC = wrapToPi(randn(pars.N, 1)*1.3);
IC = pi/2 * ones(pars.N, 1);

pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
odeoptions = odeset('RelTol', 1.0e-8,'AbsTol', 1.0e-8);
optimopts = optimoptions('fsolve', 'Display','off', 'Algorithm', 'Levenberg-Marquardt');

%% 0. Perform a full scale simulation of a FULLY CONNECTED network:
% The simple DOPRI integration: find initial conditions with fsolve
fdpars = make_fixeddegreeparameters(pars, pars.N - 1);
[t, thetas] = DOPRI_threshold(@thetaneurons, tnow, tend, IC, h, pars);
z = orderparameter(thetas);

% [t, thetas] = ode45(@(t,x) thetaneurons(t,x,pars.e,pars.K/pars.N,pars.a_n), [tnow, tend], IC, odeoptions);
% thetas = wrapToPi(thetas)';
% z = orderparameter(thetas);
disp('Small scale test done')

% The full scale simulation using the adjacency matrix:
[tfull, thetasfull] = DOPRI_simulatenetwork(tnow,tend,IC,h,fdpars);
zfull = orderparameter(thetasfull);
disp('Full scale test done')

% The mean field theory for fixed degree networks:
[T, Z] = ode45(@(t,x) MFR2(t,x,pars), [tnow, tend], gather(zfull(1)), odeoptions);
disp('Mean field test done')

% The OA mean field theory:
fdpars = prepareOAparameters(fdpars);
[TOA, ZOA, b] = OA_simulatenetwork(tnow, tend, IC, fdpars, odeoptions);
disp('OA mean field test done')

%% Plotting the results:
f_fullyconnected = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(t, abs(z), '-', 'LineWidth', 5, 'Color', '#EDB120');
plot(tfull, abs(zfull), '-', 'LineWidth', 4, 'Color', '#0072BD');
plot(T, abs(Z), '-', 'LineWidth', 3, 'Color', '#D95319');
plot(TOA, abs(ZOA), '-', 'LineWidth', 2, 'Color', '#000000');
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Fully connected network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fdpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{simple}$$', '$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest', 'Orientation','horizontal')
exportpdf(f_fullyconnected, '../Figures/InspectMeanFieldFullyConnected.pdf', export);
close(f_fullyconnected)

disp('Made fully connected network figure')

%% 1. Perform a full scale simulation of a fixed degree network:
netdegree = round(pars.N*0.3);

% The full scale simulation using the adjacency matrix:
fdpars = make_fixeddegreeparameters(pars, netdegree);
[tfull, thetasfull] = DOPRI_simulatenetwork(tnow,tend,IC,h,fdpars);
zfull = orderparameter(thetasfull);
disp('Full scale test done')

% The mean field theory for fixed degree networks:
[T, Z] = ode45(@(t,x) MFR(t,x,pars), [tnow, tend], gather(zfull(1)), odeoptions);
disp('Mean field test done')

% The OA mean field theory:
fdpars = prepareOAparameters(fdpars);
[TOA, ZOA] = OA_simulatenetwork(tnow, tend, IC, fdpars, odeoptions);
disp('OA mean field test done')

%% Plotting the results:
f_fixeddegree = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(tfull, abs(zfull), '-', 'LineWidth', 4, 'Color', '#0072BD');
plot(T, abs(Z), '-', 'LineWidth', 3, 'Color', '#D95319');
plot(TOA, abs(ZOA), '-', 'LineWidth', 2, 'Color', '#000000');
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Fixed degree network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fdpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')

legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest', 'Orientation','horizontal')
exportpdf(f_fixeddegree, '../Figures/InspectMeanFieldFixedDegree.pdf', export);
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

%% Plotting the results:
f_random = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(tfull, abs(zfull), '-', 'LineWidth', 3, 'Color', '#0072BD');
plot(TOA, abs(ZOA), '-', 'LineWidth', 2, 'Color', '#000000');
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Random network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$p$$ = %0.1f', pars.N, rdpars.meandegree, rdpars.netp), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest', 'Orientation','horizontal')
exportpdf(f_random, '../Figures/InspectMeanFieldRandom.pdf', export);
close(f_random)

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

%% Plotting the results:
f_scalefree = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(tfull, abs(zfull), '-', 'LineWidth', 3, 'Color', '#0072BD');
plot(TOA, abs(ZOA), '-k', 'LineWidth', 2, 'Color', '#000000');
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Scale-free network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$\\gamma$$ = %0.1f', pars.N, sfpars.meandegree, sfpars.degree), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest', 'Orientation','horizontal')
exportpdf(f_scalefree, '../Figures/InspectMeanFieldScaleFree.pdf', export);
close(f_scalefree)

disp('Made scale-free network figure')

% Testing the OAIC 
%% From z to ZOA
clc
IC = 0.2*randn(pars.N, 1);
IC = linspace(-2*pi, 2*pi, pars.N + 1);
IC = IC(1:end-1);
% IC = rand(pars.N, 1) * 2*pi - pi;
z = orderparameter(IC)

sfpars = make_scalefreeparameters(pars, 4);
sfpars = prepareOAparameters(sfpars);

OAIC = rand(1,sfpars.l);
for i = 1:sfpars.l
    meanthetaperdegree = IC(sfpars.degrees_in == sfpars.k(i));
%     OAIC(i) = orderparameter(meanthetaperdegree) / sfpars.P(sfpars.k(i)) * numel(meanthetaperdegree);
    OAIC(i) = sum(exp(1i*meanthetaperdegree)) / sfpars.P(sfpars.k(i));
%     OAIC(i) = mean(meanthetaperdegree);
end
Z = OAIC*sfpars.P(sfpars.k)/sfpars.N

histogram(abs(OAIC), linspace(-2*pi, 2*pi, 10), 'Normalization' , 'pdf')

% Z = orderparameter(OAIC)
% Z = ones(1,sfpars.N)*exp(1i*IC)/sfpars.N

%% From ZOA to z
clc
sfpars = make_scalefreeparameters(pars, 2.1);
sfpars = prepareOAparameters(sfpars);

OAIC = rand(1, sfpars.l);
testOAIC = orderparameter(OAIC)

Z = exp(1i*OAIC)*sfpars.P(sfpars.k)/sfpars.N

IC = zeros(pars.N, 1);
for j = 1:sfpars.l
    idx = sfpars.degrees_in == sfpars.k(j);
    IC(idx) = OAIC(j)*sfpars.P(sfpars.k(j));
end
OAIC(1)
IC(1)
testIC = orderparameter(IC)

