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
tnow = 0; tend = 3;
h = 0.01;

pars.N = 500;
pars.a_n = 0.666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
seed = 1; rng(seed);
IC = randn(pars.N, 1) + 1;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

%% 0. Perform a full scale simulation of a FULLY CONNECTED network:
[~, thetas] = DOPRI_threshold(@thetaneurons, tnow, tend, IC, h, pars);
z = orderparameter(thetas);

fixeddegreepars = make_fixeddegreeparameters(pars, pars.N - 1);
[t, thetas_full] = DOPRI_simulatenetwork(tnow,tend,IC,h,fixeddegreepars);
z_full = orderparameter(thetas_full);

% Predict order parameter from mean field theory for fixed degree networks:
MRFIC = z(1);
options = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
[T, Z] = ode45(@(t,x) MFR(t,x,pars), [tnow, tend], MRFIC, options);

%% Plotting the results:
f_fullyconnected = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(t, abs(z), '-', 'LineWidth', 2);
plot(t, abs(z_full), '-', 'LineWidth', 2);
plot(T, abs(Z), '-', 'LineWidth', 2);
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

legend('$$Z(t)_{simple}$$', '$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest')
removewhitspace();

print(f_fullyconnected, '../Figures/InspectMeanFieldFullyConnected.png', '-dpng', '-r300')

%% 1. Perform a full scale simulation of a FULLY CONNECTED network:

