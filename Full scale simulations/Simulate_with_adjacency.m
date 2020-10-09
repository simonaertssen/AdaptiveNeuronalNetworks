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
h = 0.02;

pars.N = 5000;
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

%%
% The OA mean field theory:
fixeddegreepars = prepareOAparameters(fixeddegreepars);
OAIC = ones(fixeddegreepars.l,1)*MFIC + 0.001*randn(fixeddegreepars.l,1);
[Toa, b_i] = ode45(@(t,x) MFROA(t,x,fixeddegreepars), [tnow, tend], OAIC, options);
% [timepoints, ks] = size(b_d);
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

legend('$$Z(t)_{simple}$$', '$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{OA}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest')
removewhitspace();

print(f_fullyconnected, '../Figures/InspectMeanFieldFullyConnected.png', '-dpng', '-r300')


