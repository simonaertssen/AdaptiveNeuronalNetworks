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

%% Theta model parameters:
tnow = 0; tend = 5;
h = 0.005;

pars.N = 15000;
pars.a_n = 0.666666666666666666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;

seed = 2; rng(seed);
IC = wrapToPi(randn(pars.N, 1)*1.3);
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
odeoptions = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);

%% Make a GPU init handle:
if gpuDeviceCount > 0
    d = gpuDevice(1);
    disp(d)
end
initarray = make_GPUhandle();

%% 0. Perform a full scale simulation of a FULLY CONNECTED network:
tic;
[t, thetas] = DOPRI_threshold(@thetaneurons, tnow, tend, IC, h, pars);
z = orderparameter(thetas);
toc

% tic;
% [t_ode45, theta_ode45] = ode45(@(t,x) thetaneurons(t,x,pars.e,pars.K/pars.N,pars.a_n), [tnow, tend], IC, odeoptions);
% theta_ode45 = wrapToPi(theta_ode45)';
% zode45 = orderparameter(theta_ode45)';
% toc

fdpars = make_fixeddegreeparameters(pars, pars.N - 1);
tic;
[t_full, thetas_full, A] = DOPRI_simulatenetwork(tnow, tend, IC, h, fdpars, fdpars.K);
z_full = orderparameter(thetas_full);
toc

% tic;
% [t_fullode45, thetas_fullode45] = ode45(@(t,x,K) thetaneurons_full(t,x,fdpars.K,A,fdpars.e,1/fdpars.meandegree,fdpars.a_n), [tnow, tend], IC, odeoptions);
% thetas_fullode45 = wrapToPi(thetas_fullode45)';
% z_fullode45 = orderparameter(thetas_fullode45)';
% toc


%% The mean field theory for fixed degree networks:
MFIC = gather(z(1));
[T, Z] = ode45(@(t,x) MFR2(t,x,pars), [tnow, tend], MFIC, odeoptions);

fdpars = prepareOAparameters(fdpars);
OAIC = ones(fdpars.l,1)*MFIC + 0.001*randn(fdpars.l,1);
[Toa, b_i] = ode45(@(t,x) MFROA(t,x,fdpars), [tnow, tend], OAIC, odeoptions);
Zoa = gather(initarray(b_i) * fdpars.P(fdpars.k)/fdpars.N);

%% Plotting:
f_CPW = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(Toa, abs(Zoa), '-r', 'LineWidth', 3);
plot(T, abs(Z), ':k', 'LineWidth', 2);

plot(t, abs(z), '-', 'LineWidth', 4);
% plot(t_ode45, abs(zode45), '-', 'LineWidth', 3);
plot(t_full, abs(z_full), '-', 'LineWidth', 2);
% plot(t_fullode45, abs(z_fullode45), '-', 'LineWidth', 1);
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

legend('$$\overline{Z(t)}_{\rm MF}$$', '$$\overline{Z(t)}_{\rm OA}$$', '$$Z(t)_{\rm DOPRI}$$', '$$Z(t)_{A_{ij},\rm DOPRI}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest')
title(sprintf('\\bf Fully connected network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fdpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')
exportpdf(f_CPW, '../Figures/ValidateCPWcycle.pdf', true);

close(f_CPW)

disp('Made fully connected network figure')