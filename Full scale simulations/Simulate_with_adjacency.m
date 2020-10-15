clear all; close all; clc;
% Simulate a full scale fixed degree network and test whether the results
% are correct, with respect to the order parameters.

%% Setup
addpath('../Functions');
addpath('../Mean Field Reductions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 15;
labelfont = 15;
export = true;

%% Make a GPU init handle:
if gpuDeviceCount > 0
    d = gpuDevice(gpuDeviceCount-1);
    disp(d)
end
initarray = make_GPUhandle();

%% Theta model parameters:
tnow = 0; tend = 5;
h = 0.001;

pars.N = 15000;
pars.a_n = 0.666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;

seed = 1; rng(seed);
IC = wrapToPi(randn(pars.N, 1)*1.3);
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
odeoptions = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6, 'NormControl','on');

%% 0. Perform a full scale simulation of a FULLY CONNECTED network:
% The simple DOPRI integration:
[t, thetas] = DOPRI_threshold(@thetaneurons, tnow, tend, IC, h, pars);
z = orderparameter(thetas);
% [t, thetas] = ode113(@(t,x) thetaneurons(t,x,pars.e,pars.K/pars.N,pars.a_n), [tnow, tend], IC, odeoptions);
% thetas = wrapToPi(thetas)';
% z = orderparameter(thetas);
disp('Small scale test done')

% The full scale simulation using the adjacency matrix:
fdpars = make_fixeddegreeparameters(pars, pars.N - 1);
[t_full, thetas_full] = DOPRI_simulatenetwork(tnow,tend,IC,h,fdpars);
z_full = orderparameter(thetas_full);
% fdpars = make_fixeddegreeparameters(pars, pars.N - 1);
% A = initarray(adjacencymatrix(fdpars.degrees_in, fdpars.degrees_out));
% [t_full, thetas_full] = ode113(@(t,x,K) thetaneurons_full(t,x,fdpars.K,A,fdpars.e,1/fdpars.meandegree,fdpars.a_n), [tnow, tend], IC, odeoptions);
% thetas_full = wrapToPi(thetas_full)';
% z_full = orderparameter(thetas_full)';
disp('Full scale test done')

% The mean field theory for fixed degree networks:
MFIC = gather(z(1));
[T, Z] = ode45(@(t,x) MFR2(t,x,pars), [tnow, tend], MFIC, odeoptions);
disp('Mean field test done')

% The OA mean field theory:
fdpars = prepareOAparameters(fdpars);
OAIC = ones(fdpars.l,1)*MFIC;
[Toa, b_i] = ode45(@(t,x) MFROA(t,x,fdpars), [tnow, tend], OAIC, odeoptions);
Zoa = gather(initarray(b_i) * fdpars.P(fdpars.k)/fdpars.N);
disp('OA mean field test done')


%% Plotting the results:
f_fullyconnected = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(t, abs(z), '-', 'LineWidth', 5, 'Color', '#EDB120');
plot(t_full, abs(z_full), '-', 'LineWidth', 4, 'Color', '#0072BD');
plot(T, abs(Z), '-', 'LineWidth', 3, 'Color', '#D95319');
plot(Toa, abs(Zoa), '-', 'LineWidth', 2, 'Color', '#000000');
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Fully connected network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fdpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{simple}$$', '$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest', 'Orientation','horizontal')
exportpdf(f_fullyconnected, '../Figures/InspectMeanFieldFullyConnected.pdf', export);
close(f_fullyconnected)

disp('Made fully connected network figure')

%% 1. Perform a full scale simulation of a fixed degree network:
% The full scale simulation using the adjacency matrix:
netdegree = round(pars.N*0.3);
fdpars = make_fixeddegreeparameters(pars, netdegree);
[t_full, thetas_full] = DOPRI_simulatenetwork(tnow,tend,IC,h,fdpars);
z_full = orderparameter(thetas_full);
% A = initarray(adjacencymatrix(fdpars.degrees_in, fdpars.degrees_out));
% [t_full, thetas_full] = ode113(@(t,x,K) thetaneurons_full(t,x,fdpars.K,A,fdpars.e,1/fdpars.meandegree,fdpars.a_n), [tnow, tend], IC, odeoptions);
% thetas_full = wrapToPi(thetas_full)';
% z_full = orderparameter(thetas_full)';
disp('Full scale test done')

% The mean field theory for fixed degree networks:
MFIC = gather(z_full(1));
[T, Z] = ode45(@(t,x) MFR(t,x,pars), [tnow, tend], MFIC, odeoptions);
disp('Mean field test done')

% The OA mean field theory:
fdpars = prepareOAparameters(fdpars);
OAIC = ones(fdpars.l,1)*MFIC;
[Toa, b_i] = ode45(@(t,x) MFROA(t,x,fdpars), [tnow, tend], OAIC, odeoptions);
Zoa = gather(initarray(b_i) * fdpars.P(fdpars.k)/fdpars.N);
disp('OA mean field test done')

%% Plotting the results:
f_fixeddegree = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(t_full, abs(z_full), '-', 'LineWidth', 4, 'Color', '#0072BD');
plot(T, abs(Z), '-', 'LineWidth', 3, 'Color', '#D95319');
plot(Toa, abs(Zoa), '-', 'LineWidth', 2, 'Color', '#000000');
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Fixed degree network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fdpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')

legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest', 'Orientation','horizontal')
exportpdf(f_fixeddegree, '../Figures/InspectMeanFieldFixedDegree.pdf', export);
close(f_fixeddegree)

disp('Made fixed degree network figure')

%% 2. Perform a full scale simulation of a random network:
% The full scale simulation using the adjacency matrix:
netp = 0.3;
rdpars = make_randomparameters(pars, netp);
[t_full, thetas_full] = DOPRI_simulatenetwork(tnow,tend,IC,h,rdpars);
z_full = orderparameter(thetas_full);
% A = initarray(adjacencymatrix(rdpars.degrees_in, rdpars.degrees_out));
% [t_full, thetas_full] = ode113(@(t,x,K) thetaneurons_full(t,x,rdpars.K,A,rdpars.e,1/rdpars.meandegree,rdpars.a_n), [tnow, tend], IC, odeoptions);
% thetas_full = wrapToPi(thetas_full)';
% z_full = orderparameter(thetas_full)';
disp('Full scale test done')

% The OA mean field theory:
rdpars = prepareOAparameters(rdpars);
OAIC = ones(rdpars.l,1)*gather(z_full(1));
[Toa, b_i] = ode45(@(t,x) MFROA(t,x,rdpars), [tnow, tend], OAIC, odeoptions);
Zoa = gather(initarray(b_i) * rdpars.P(rdpars.k)/rdpars.N);
disp('OA mean field test done')

%% Plotting the results:
f_random = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(t_full, abs(z_full), '-', 'LineWidth', 3, 'Color', '#0072BD');
plot(Toa, abs(Zoa), '-', 'LineWidth', 2, 'Color', '#000000');
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Random network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$p$$ = %0.1f', pars.N, rdpars.meandegree, rdpars.netp), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest', 'Orientation','horizontal')
exportpdf(f_random, '../Figures/InspectMeanFieldRandom.pdf', export);
close(f_random)

disp('Made random network figure')

%% 3. Perform a full scale simulation of a scale-free network:
% The full scale simulation using the adjacency matrix:
IC = wrapToPi(randn(pars.N, 1)*0.5);
degree = 4;
sfpars = make_scalefreeparameters(pars, degree);
[t_full, thetas_full] = DOPRI_simulatenetwork(tnow,tend,IC,h,sfpars);
z_full = orderparameter(thetas_full);
% A = initarray(adjacencymatrix(sfpars.degrees_in, sfpars.degrees_out));
% [t_full, thetas_full] = ode113(@(t,x,K) thetaneurons_full(t,x,sfpars.K,A,sfpars.e,1/sfpars.meandegree,sfpars.a_n), [tnow, tend], IC, odeoptions);
% thetas_full = wrapToPi(thetas_full)';
% z_full = orderparameter(thetas_full)';
disp('Full scale test done')

% The OA mean field theory:
sfpars = prepareOAparameters(sfpars);
OAIC = ones(sfpars.l,1)*gather(z_full(1));
[Toa, b_i] = ode45(@(t,x) MFROA(t,x,sfpars), [tnow, tend], OAIC, odeoptions);
Zoa = gather(initarray(b_i) * sfpars.P(sfpars.k)/sfpars.N);
disp('OA mean field test done')

%% Plotting the results:
f_scalefree = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(t_full, abs(z_full), '-', 'LineWidth', 3, 'Color', '#0072BD');
plot(Toa, abs(Zoa), '-k', 'LineWidth', 2, 'Color', '#000000');
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Scale-free network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$\\gamma$$ = %0.1f', pars.N, sfpars.meandegree, sfpars.degree), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest', 'Orientation','horizontal')
exportpdf(f_scalefree, '../Figures/InspectMeanFieldScaleFree.pdf', export);
close(f_scalefree)

disp('Made scale-free network figure')
