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

%% Theta model parameters:
tnow = 0; tend = 5;
h = 0.001;

pars.N = 1000;
pars.a_n = 0.666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;

seed = 1; rng(seed);
IC = wrapToPi(randn(pars.N, 1)*0.8);
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

%% 0. Perform a full scale simulation of a FULLY CONNECTED network:
% The simple DOPRI integration:
tic;
[t, theta] = DOPRI_threshold(@thetaneurons, tnow, tend, IC, h, pars);
z = orderparameter(theta);
toc

tic;
options = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6, 'NormControl','on');
[tode45, theta_ode45] = ode113(@(t,x) thetaneurons(t,x,pars.e,pars.K/pars.N,pars.a_n), [tnow, tend], IC, options);
theta_ode45 = wrapToPi(theta_ode45);
zode45 = orderparameter(theta_ode45');
toc

%% The ode45 implementation
% t = []; thetas = [];
% tstart = tnow; tfinal = tend; loopcrit = true;
% options = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6, 'Events', @theta_neuron_threshold);
% 
% while loopcrit
%     [t,x] = ode45(@(t,x) thetaneurons(t,x,pars.e,pars.K/pars.N,pars.a_n), [tstart, tfinal], IC, options);
%     t = [t, t];
%     thetas = [thetas, x];
%     
%     tstart = t(end);
%     if tstart >= tend
%         loopcrit = false;
%     else
%         tspan = [tstart; tfinal];
%         IC = -pi;
%     end
% end
% z = orderparameter(thetas);

%% The mean field theory for fixed degree networks:
MFIC = gather(z(1));
options = odeset('RelTol', 1.0e-8,'AbsTol', 1.0e-8);
[T, Z] = ode45(@(t,x) MFR2(t,x,pars), [tnow, tend], MFIC, options);

%% Plotting:
f_fullyconnected = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;

xlim([tnow, tend]); ylim([0, 1])
plot(t, abs(z), '-', 'LineWidth', 2);
plot(tode45, abs(zode45), '-', 'LineWidth', 2);
plot(T, abs(Z), '-', 'LineWidth', 2);
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\vert Z (t) \vert$','Interpreter','latex', 'FontSize', labelfont)

legend('$$Z(t)_{\rm DOPRI}$$', '$$Z(t)_{\rm ode45}$$', '$$\overline{Z(t)}_{\rm MF}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southwest')
removewhitspace();

% close(f_fullyconnected)
disp('Made fully connected network figure')

%% Functions:
function [value,isterminal,direction] = theta_neuron_threshold(t,x)
    value = t - 10*pi; % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end