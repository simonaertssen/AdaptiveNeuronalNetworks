% In this script we will be studying the reponse of the theta neuron model
% on different types of input current.

clear all; close all; clc;

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

tnow = 0; tend = 10;
F = @thetaneuron; h = 0.001;

fexcite = figure('Position', [50 800 1000 300]);
titlefont = 15;
labelfont = 15;
m = 2; n = 3;

%% Excitability:
% The theta neuron model is of type 1 excitability, meaning that the frequency 
% of spikes goes up as the input current increases.
I = @excitabilitycurrent;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, -pi, h, I);
drawthetas = spikesNaN(thetas);

imrow(1) = subplot(m,n,1); hold on; box on;
removewhitspace();

title("Class 1 excitability", 'FontSize', titlefont, 'FontName', 'SansSerif');
yyaxis left
ylim([-pi - 1.5, pi + 0.3]);
plot(t, thetas, ':k', 'LineWidth', 1)
plot(t, drawthetas, '-', 'LineWidth', 2.5, 'color', '#0072BD');
ylabel('$\theta_i$','Interpreter','latex', 'FontSize', labelfont);

yyaxis right
maxy = max(I(t));
plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';
set(gca,'YTick', 0:maxy:maxy, 'YLim', [-250, maxy*8]);

removewhitspace();


%% Spiking behaviour:
% The neuron spikes when the input current > 0
I = @spikecurrent;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, -pi, h, I);
drawthetas = spikesNaN(thetas);

imrow(2) = subplot(m,n,2); hold on; box on;
removewhitspace();

title("Spiking", 'FontSize', titlefont, 'FontName', 'SansSerif');
yyaxis left
ylim([-pi - 1.5, pi + 0.3]);
plot(t, thetas, ':k', 'LineWidth', 1)
plot(t, drawthetas, '-', 'LineWidth', 2.5, 'color', '#0072BD');

yyaxis right
maxy = max(I(t));
plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';
set(gca,'YTick', 0:maxy:maxy, 'YLim', [-1, maxy*8]);

removewhitspace();


%% Bursting behaviour:
% The neuron spikes when the input current > 0
I = @burstcurrent;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, -pi, h, I);
drawthetas = spikesNaN(thetas);

imrow(3) = subplot(m,n,3); hold on; box on;
removewhitspace();

title("Bursting", 'FontSize', titlefont, 'FontName', 'SansSerif');
yyaxis left
ylim([-pi - 1.5, pi + 0.3]);
plot(t, thetas, ':k', 'LineWidth', 1)
plot(t, drawthetas, '-', 'LineWidth', 2.5, 'color', '#0072BD');

yyaxis right
maxy = max(I(t));
plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
ylim([-1, max(I(t))*5]);
ylabel('$I$','Interpreter','latex', 'FontSize', labelfont)
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';
set(gca,'YTick', 0:maxy:maxy, 'YLim', [-200, maxy*8])

removewhitspace();


%% Rebound behaviour:
% The neuron spikes when the input current > 0
I = @rebound;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, -pi, h, I);
drawthetas = spikesNaN(thetas);

imrow(3) = subplot(m,n,3); hold on; box on;
removewhitspace();

title("Bursting", 'FontSize', titlefont, 'FontName', 'SansSerif');
yyaxis left
ylim([-pi - 1.5, pi + 0.3]);
plot(t, thetas, ':k', 'LineWidth', 1)
plot(t, drawthetas, '-', 'LineWidth', 2.5, 'color', '#0072BD');

yyaxis right
maxy = max(I(t));
plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
ylim([-1, max(I(t))*5]);
ylabel('$I$','Interpreter','latex', 'FontSize', labelfont)
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';
set(gca,'YTick', 0:maxy:maxy, 'YLim', [-200, maxy*8])

removewhitspace();


%% Functions:
function I = excitabilitycurrent(t)
    I = max(t/10, power(t,3));
end

function I = spikecurrent(t)
    I = zeros(size(t));
    spike = 5;
    I(t > 1 & t < 2) = spike;
    I(t > 3 & t < 4) = spike;
    
    I(t > 5 & t < 6) = spike;
    I(t > 7 & t < 8) = spike;
end

function I = burstcurrent(t)
    I = zeros(size(t));
    burst = 500;
    I(t > 2 & t < 3) = burst;
    I(t > 6 & t < 7) = burst;
end

function I = burstcurrent(t)
    I = zeros(size(t));
    burst = 500;
    I(t > 2 & t < 3) = burst;
    I(t > 6 & t < 7) = burst;
end
