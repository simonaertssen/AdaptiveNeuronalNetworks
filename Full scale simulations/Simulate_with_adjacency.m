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

%% Make a GPU init handle:
if gpuDeviceCount > 0
    d = gpuDevice(gpuDeviceCount-1);
    disp(d)
end
initarray = make_GPUhandle();

%% Theta model parameters:
tnow = 0; tend = 10;
h = 0.001;

pars.N = 10000;
pars.a_n = 0.666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
% pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;

seed = 1; rng(seed);
IC = randn(pars.N, 1) + 1;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
odeoptions = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6, 'NormControl','on');

%% 0. Perform a full scale simulation of a FULLY CONNECTED network:
% The simple DOPRI integration:
% [~, thetas] = DOPRI_threshold(@thetaneurons, tnow, tend, IC, h, pars);
% z = orderparameter(thetas);
[t, thetas] = ode113(@(t,x) thetaneurons(t,x,pars.e,pars.K/pars.N,pars.a_n), [tnow, tend], IC, odeoptions);
thetas = wrapToPi(thetas);
z = orderparameter(thetas');

% The full scale simulation using the adjacency matrix:
% fixeddegreepars = make_fixeddegreeparameters(pars, pars.N - 1);
% [t, thetas_full] = DOPRI_simulatenetwork(tnow,tend,IC,h,fixeddegreepars);
% z_full = orderparameter(thetas_full);
fdpars = make_fixeddegreeparameters(pars, pars.N - 1);
A = initarray(adjacencymatrix(fdpars.degrees_in, fdpars.degrees_out));
[t_full, thetas_full] = ode113(@(t,x,K) thetaneurons_full(t,x,fdpars.K,A,fdpars.e,1/fdpars.meandegree,fdpars.a_n), [tnow, tend], IC, odeoptions);
thetas_full = wrapToPi(thetas_full);
z_full = orderparameter(thetas_full');

% The mean field theory for fixed degree networks:
MFIC = gather(z(1));
options = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
[T, Z] = ode45(@(t,x) MFR(t,x,pars), [tnow, tend], MFIC, options);

% The OA mean field theory:
fdpars = prepareOAparameters(fdpars);
OAIC = ones(fdpars.l,1)*MFIC + 0.001*randn(fdpars.l,1);
[Toa, b_i] = ode45(@(t,x) MFROA(t,x,fdpars), [tnow, tend], OAIC, options);
Zoa = b_i * fdpars.P(fdpars.k)/fdpars.N;

%% Plotting the results:
f_fullyconnected = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(t, abs(z), '-', 'LineWidth', 2);
plot(t_full, abs(z_full), '-', 'LineWidth', 2);
plot(T, abs(Z), '-', 'LineWidth', 2);
plot(Toa, abs(Zoa), '-', 'LineWidth', 2);
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Fully connected network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fdpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{simple}$$', '$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{OA}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest')
removewhitspace();

print(f_fullyconnected, '../Figures/InspectMeanFieldFullyConnected.png', '-dpng', '-r300')
close(f_fullyconnected)
disp('Made fully connected network figure')

%% 1. Perform a full scale simulation of a fixed degree network:
% The full scale simulation using the adjacency matrix:
netdegree = round(pars.N*0.3);
fdpars = make_fixeddegreeparameters(pars, netdegree);
% [t, thetas_full] = DOPRI_simulatenetwork(tnow,tend,IC,h,fdpars);
% z_full = orderparameter(thetas_full);
A = initarray(adjacencymatrix(fdpars.degrees_in, fdpars.degrees_out));
[t_full, thetas_full] = ode113(@(t,x,K) thetaneurons_full(t,x,fdpars.K,A,fdpars.e,1/fdpars.meandegree,fdpars.a_n), [tnow, tend], IC, odeoptions);
thetas_full = wrapToPi(thetas_full);
z_full = orderparameter(thetas_full');

% The mean field theory for fixed degree networks:
MFIC = gather(z_full(1));
options = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
[T, Z] = ode45(@(t,x) MFR(t,x,pars), [tnow, tend], MFIC, options);

% The OA mean field theory:
fdpars = prepareOAparameters(fdpars);
OAIC = ones(fdpars.l,1)*MFIC + 0.001*randn(fdpars.l,1);
[Toa, b_i] = ode45(@(t,x) MFROA(t,x,fdpars), [tnow, tend], OAIC, options);
Zoa = b_i * fdpars.P(fdpars.k)/fdpars.N;

%% Plotting the results:
f_fixeddegree = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(t_full, abs(z_full), '-', 'LineWidth', 2);
plot(T, abs(Z), '-', 'LineWidth', 2);
plot(Toa, abs(Zoa), '-', 'LineWidth', 2);
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Fixed degree network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fdpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')

legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{OA}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest')
removewhitspace();

print(f_fixeddegree, '../Figures/InspectMeanFieldFixedDegree.png', '-dpng', '-r300')
close(f_fixeddegree)
disp('Made fixed degree network figure')

%% 2. Perform a full scale simulation of a random network:
% The full scale simulation using the adjacency matrix:
netp = 0.3;
rdpars = make_randomparameters(pars, netp);
% [t, thetas_full] = DOPRI_simulatenetwork(tnow,tend,IC,h,randompars);
% z_full = orderparameter(thetas_full);
A = initarray(adjacencymatrix(rdpars.degrees_in, rdpars.degrees_out));
[t_full, thetas_full] = ode113(@(t,x,K) thetaneurons_full(t,x,rdpars.K,A,rdpars.e,1/rdpars.meandegree,rdpars.a_n), [tnow, tend], IC, odeoptions);
thetas_full = wrapToPi(thetas_full);
z_full = orderparameter(thetas_full');

% The mean field theory for fixed degree networks:
MFIC = gather(z_full(1));
options = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
[T, Z] = ode45(@(t,x) MFR(t,x,pars), [tnow, tend], MFIC, options);

% The OA mean field theory:
rdpars = prepareOAparameters(rdpars);
OAIC = ones(rdpars.l,1)*MFIC + 0.001*randn(rdpars.l,1);
[Toa, b_i] = ode45(@(t,x) MFROA(t,x,rdpars), [tnow, tend], OAIC, options);
Zoa = b_i * rdpars.P(rdpars.k)/rdpars.N;

%% Plotting the results:
f_random = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(t_full, abs(z_full), '-', 'LineWidth', 2);
plot(T, abs(Z), '-', 'LineWidth', 2);
plot(Toa, abs(Zoa), '-', 'LineWidth', 2);
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Random network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$p$$ = %0.1f', pars.N, rdpars.meandegree, rdpars.netp), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{OA}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest')
removewhitspace();

print(f_random, '../Figures/InspectMeanFieldRandom.png', '-dpng', '-r300')
close(f_random)
disp('Made random network figure')

%% 3. Perform a full scale simulation of a scale-free network:
% The full scale simulation using the adjacency matrix:
degree = 3;
sfpars = make_scalefreeparameters(pars, degree);
% [t, thetas_full] = DOPRI_simulatenetwork(tnow,tend,IC,h,scalefreepars);
% z_full = orderparameter(thetas_full);
A = initarray(adjacencymatrix(sfpars.degrees_in, sfpars.degrees_out));
[t_full, thetas_full] = ode113(@(t,x,K) thetaneurons_full(t,x,sfpars.K,A,sfpars.e,1/sfpars.meandegree,sfpars.a_n), [tnow, tend], IC, odeoptions);
thetas_full = wrapToPi(thetas_full);
z_full = orderparameter(thetas_full');

% The mean field theory for fixed degree networks:
MFIC = gather(z_full(1));
options = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
[T, Z] = ode45(@(t,x) MFR(t,x,pars), [tnow, tend], MFIC, options);

% The OA mean field theory:
sfpars = prepareOAparameters(sfpars);
OAIC = ones(sfpars.l,1)*MFIC + 0.001*randn(sfpars.l,1);
[Toa, b_i] = ode45(@(t,x) MFROA(t,x,sfpars), [tnow, tend], OAIC, options);
Zoa = b_i * sfpars.P(sfpars.k)/sfpars.N;

%% Plotting the results:
f_scalefree = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(t_full, abs(z_full), '-', 'LineWidth', 2);
plot(T, abs(Z), '-', 'LineWidth', 2);
plot(Toa, abs(Zoa), '-', 'LineWidth', 2);
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Scale-free network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$\\gamma$$ = %0.1f', pars.N, sfpars.meandegree, sfpars.degree), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{OA}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest')
removewhitspace();

print(f_scalefree, '../Figures/InspectMeanFieldScaleFree.png', '-dpng', '-r300')
close(f_scalefree)
disp('Made scale-free network figure')
