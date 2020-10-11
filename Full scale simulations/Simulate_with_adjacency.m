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
    d = gpuDevice(gpuDeviceCount);
    disp(d)
end
initarray = make_GPUhandle();

%% Theta model parameters:
tnow = 0; tend = 2;
h = 0.005;

pars.N = 1000;
pars.a_n = 0.666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
seed = 1; rng(seed);
IC = randn(pars.N, 1) + 1;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

%% 0. Perform a full scale simulation of a FULLY CONNECTED network:
% The simple DOPRI integration:
[~, thetas] = DOPRI_threshold(@thetaneurons, tnow, tend, IC, h, pars);
z = orderparameter(thetas);

% The full scale simulation using the adjacency matrix:
fixeddegreepars = make_fixeddegreeparameters(pars, pars.N - 1);
[t, thetas_full] = DOPRI_simulatenetwork(tnow,tend,IC,h,fixeddegreepars);
z_full = orderparameter(thetas_full);

% The mean field theory for fixed degree networks:
MFIC = z(1);
options = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
[T, Z] = ode45(@(t,x) MFR(t,x,pars), [tnow, tend], MFIC, options);

% The OA mean field theory:
fixeddegreepars = prepareOAparameters(fixeddegreepars);
OAIC = ones(fixeddegreepars.l,1)*MFIC + 0.001*randn(fixeddegreepars.l,1);
[Toa, b_i] = ode45(@(t,x) MFROA(t,x,fixeddegreepars), [tnow, tend], OAIC, options);
Zoa = b_i * fixeddegreepars.P(fixeddegreepars.k)/fixeddegreepars.N;

%% Plotting the results:
f_fullyconnected = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(t, abs(z), '-', 'LineWidth', 2);
plot(t, abs(z_full), '-', 'LineWidth', 2);
plot(T, abs(Z), '-', 'LineWidth', 2);
plot(Toa, abs(Zoa), '-', 'LineWidth', 2);
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Fully connected network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fixeddegreepars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{simple}$$', '$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{OA}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest')
removewhitspace();

print(f_fullyconnected, '../Figures/InspectMeanFieldFullyConnected.png', '-dpng', '-r300')
close(f_fullyconnected)

%% 1. Perform a full scale simulation of a fixed degree network:
% The full scale simulation using the adjacency matrix:
netdegree = round(pars.N*0.3);
fixeddegreepars = make_fixeddegreeparameters(pars, netdegree);
[t, thetas_full] = DOPRI_simulatenetwork(tnow,tend,IC,h,fixeddegreepars);
z_full = orderparameter(thetas_full);

% The mean field theory for fixed degree networks:
MFIC = z_full(1);
options = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
[T, Z] = ode45(@(t,x) MFR(t,x,pars), [tnow, tend], MFIC, options);

% The OA mean field theory:
fixeddegreepars = prepareOAparameters(fixeddegreepars);
OAIC = ones(fixeddegreepars.l,1)*MFIC + 0.001*randn(fixeddegreepars.l,1);
[Toa, b_i] = ode45(@(t,x) MFROA(t,x,fixeddegreepars), [tnow, tend], OAIC, options);
Zoa = b_i * fixeddegreepars.P(fixeddegreepars.k)/fixeddegreepars.N;

%% Plotting the results:
f_fixeddegree = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(t, abs(z_full), '-', 'LineWidth', 2);
plot(T, abs(Z), '-', 'LineWidth', 2);
plot(Toa, abs(Zoa), '-', 'LineWidth', 2);
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Fixed degree network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fixeddegreepars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')

legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{OA}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest')
removewhitspace();

print(f_fixeddegree, '../Figures/InspectMeanFieldFixedDegree.png', '-dpng', '-r300')
close(f_fixeddegree)

%% 2. Perform a full scale simulation of a random network:
% The full scale simulation using the adjacency matrix:
netp = 0.3;
randompars = make_randomparameters(pars, netp);
[t, thetas_full] = DOPRI_simulatenetwork(tnow,tend,IC,h,randompars);
z_full = orderparameter(thetas_full);

% The mean field theory for fixed degree networks:
MFIC = z_full(1);
options = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
[T, Z] = ode45(@(t,x) MFR(t,x,pars), [tnow, tend], MFIC, options);

% The OA mean field theory:
randompars = prepareOAparameters(randompars);
OAIC = ones(randompars.l,1)*MFIC + 0.001*randn(randompars.l,1);
[Toa, b_i] = ode45(@(t,x) MFROA(t,x,randompars), [tnow, tend], OAIC, options);
Zoa = b_i * randompars.P(randompars.k)/randompars.N;

%% Plotting the results:
f_random = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(t, abs(z_full), '-', 'LineWidth', 2);
plot(T, abs(Z), '-', 'LineWidth', 2);
plot(Toa, abs(Zoa), '-', 'LineWidth', 2);
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Random network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$p$$ = %0.1f', pars.N, randompars.meandegree, randompars.netp), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{OA}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest')
removewhitspace();

print(f_random, '../Figures/InspectMeanFieldRandom.png', '-dpng', '-r300')
close(f_random)

%% 3. Perform a full scale simulation of a scale-free network:
% The full scale simulation using the adjacency matrix:
degree = 3;
scalefreepars = make_scalefreeparameters(pars, degree);
[t, thetas_full] = DOPRI_simulatenetwork(tnow,tend,IC,h,scalefreepars);
z_full = orderparameter(thetas_full);

% The mean field theory for fixed degree networks:
MFIC = z_full(1);
options = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
[T, Z] = ode45(@(t,x) MFR(t,x,pars), [tnow, tend], MFIC, options);

% The OA mean field theory:
scalefreepars = prepareOAparameters(scalefreepars);
OAIC = ones(scalefreepars.l,1)*MFIC + 0.001*randn(scalefreepars.l,1);
[Toa, b_i] = ode45(@(t,x) MFROA(t,x,scalefreepars), [tnow, tend], OAIC, options);
Zoa = b_i * scalefreepars.P(scalefreepars.k)/scalefreepars.N;

%% Plotting the results:
f_scalefree = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(t, abs(z_full), '-', 'LineWidth', 2);
plot(T, abs(Z), '-', 'LineWidth', 2);
plot(Toa, abs(Zoa), '-', 'LineWidth', 2);
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

title(sprintf('\\bf Scale-free network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$\\gamma$$ = %0.1f', pars.N, scalefreepars.meandegree, scalefreepars.degree), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{OA}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest')
removewhitspace();

print(f_scalefree, '../Figures/InspectMeanFieldScaleFree.png', '-dpng', '-r300')
close(f_scalefree)
