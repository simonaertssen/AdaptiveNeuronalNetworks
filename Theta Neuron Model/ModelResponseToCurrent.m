% In this script we will be studying the reponse of the theta neuron model
% on different types of input current.

clear all; close all; clc;

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

tnow = 0; tend = 9;
F = @thetaneuron; h = 0.001;

%% Excitability:
% The theta neuron model is of type 1 excitability, meaning that the frequency 
% of spikes goes up as the input current increases.
I = @excitabilitycurrent;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, -pi, h, I);
drawthetas = spikesNaN(thetas);

fexcite = figure('Position', [50 800 1000 200]); hold on; box on
yyaxis left
ylim([-pi - 0.1, pi + 0.1]);
plot(t, thetas, ':k', 'LineWidth', 1)
plot(t, drawthetas, '-', 'LineWidth', 2.5, 'color', '#0072BD');
ylabel('$\theta_i$','Interpreter','latex', 'FontSize', 20)

yyaxis right
plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
ylim([-100, max(I(t))*5]);
ylabel('$I$','Interpreter','latex', 'FontSize', 20)
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';
set(gca, 'xticklabel', []);

removewhitspace();


%% Spiking behaviour:
% The neuron spikes when the input current > 0
I = @spikecurrent;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, -pi, h, I);
drawthetas = spikesNaN(thetas);

fspiking = figure('Position', [50 800 1000 200]); hold on; box on
yyaxis left
ylim([-pi - 0.1, pi + 0.1]);
plot(t, thetas, ':k', 'LineWidth', 1)
plot(t, drawthetas, '-', 'LineWidth', 2.5, 'color', '#0072BD');
ylabel('$\theta_i$','Interpreter','latex', 'FontSize', 20)

yyaxis right
plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
ylim([-1, max(I(t))*5]);
ylabel('$I$','Interpreter','latex', 'FontSize', 20)
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';

removewhitspace();


%% Bursting behaviour:
% The neuron spikes when the input current > 0
I = @burstcurrent;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, -pi, h, I);
drawthetas = spikesNaN(thetas);

fexcite = figure('Position', [50 800 1000 200]); hold on; box on
yyaxis left
ylim([-pi - 0.1, pi + 0.1]);
plot(t, thetas, ':k', 'LineWidth', 1)
plot(t, drawthetas, '-', 'LineWidth', 2.5, 'color', '#0072BD');
ylabel('$\theta_i$','Interpreter','latex', 'FontSize', 20)

yyaxis right
plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
ylim([-1, max(I(t))*5]);
ylabel('$I$','Interpreter','latex', 'FontSize', 20)
xlabel('$t$','Interpreter','latex', 'FontSize', 20)
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];


%% Functions:
function I = excitabilitycurrent(t)
    I = max(t/10, power(t,2.5));
end

function I = spikecurrent(t)
    I = zeros(size(t));
    spike = 50;
    I(t > 1 & t < 2) = spike;
    I(t > 3 & t < 4) = spike;
    
    burst = 500;
    I(t > 5 & t < 6) = burst;
    I(t > 7 & t < 8) = burst;
end

function I = burstcurrent(t)
    I = zeros(size(t));
    bump = -500;
    I(t > 2 & t < 3) = bump;
    I(t > 6 & t < 7) = bump;
end

% For a tight layout:
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];

