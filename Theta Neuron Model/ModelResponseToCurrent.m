% In this script we will be studying the reponse of the theta neuron model
% on different types of input current.

clear all; close all; clc;

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

tnow = 0; tend = 10;
F = @thetaneuron; h = 0.001;

fexcite = figure('Renderer', 'painters', 'Position', [50 800 1000 200]);
titlefont = 15;
labelfont = 15;
m = 1; n = 3;

%% Excitability:
% The theta neuron model is of type 1 excitability, meaning that the frequency 
% of spikes goes up as the input current increases.
I = @excitabilitycurrent;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, -pi, h, I);
NaNthetas = spikesNaN(thetas);
drawthetas = 1 + cos(thetas);

imrow(1) = subplot(m,n,1); hold on; box on;

title("Class 1 excitability", 'FontSize', titlefont, 'FontName', 'SansSerif');
yyaxis left
ylim([-0.5, 2.1]);
plot(t, drawthetas, '-', 'LineWidth', 2, 'color', '#0072BD');
ylabel('$1 + \cos\theta$','Interpreter','latex', 'FontSize', labelfont);
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)

yyaxis right
maxy = max(I(t));
plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';
set(gca,'YTick', 0:maxy:maxy, 'YLim', [-10, maxy*8]);



%% Spiking behaviour:
% The neuron spikes when the input current > 0
I = @spikecurrent;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, -pi, h, I);
drawthetas = 1 + cos(thetas);

imrow(2) = subplot(m,n,2); hold on; box on;

title("Spiking", 'FontSize', titlefont, 'FontName', 'SansSerif');
yyaxis left
ylim([-0.5, 2.1]);
plot(t, drawthetas, '-', 'LineWidth', 2, 'color', '#0072BD');
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)

yyaxis right
maxy = max(I(t));
plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';
set(gca,'YTick', 0:maxy:maxy, 'YLim', [-2, maxy*8]);


%% Bursting behaviour:
% The neuron spikes when the input current > 0
I = @burstcurrent;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, -pi, h, I);
drawthetas = 1 + cos(thetas);

imrow(3) = subplot(m,n,3); hold on; box on;

title("Bursting", 'FontSize', titlefont, 'FontName', 'SansSerif');
yyaxis left
ylim([-0.5, 2.1]);
plot(t, drawthetas, '-', 'LineWidth', 2, 'color', '#0072BD');
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)

yyaxis right
maxy = max(I(t));
plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
ylim([-1, max(I(t))*5]);
ylabel('$I$','Interpreter','latex', 'FontSize', labelfont)
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';
set(gca,'YTick', 0:maxy:maxy, 'YLim', [-200, maxy*8])


%% Save the figure:
print(fexcite, '../Figures/ThetaNeuronResponseToCurrent.png', '-dpng', '-r300')
%exportgraphics(fexcite, '../Figures/ThetaNeuronResponseToCurrent.png', 'Resolution', 300')


%% Investigate the frequency - current curve:
% We know that T = pi/sqrt(I)

fI = figure('Renderer', 'painters', 'Position', [50 800 500 200]); hold on; box on;

I = @fI_current;

tend = 50;
[t, thetas] = DOPRI_singleneuron(F, 0, tend, -pi, h, I);
NaNthetas = spikesNaN(thetas);
drawthetas = 1 + cos(thetas);

idx = [1, find(isnan(NaNthetas))];
NaNnum = length(idx)-1;
frequencies = diff(t(idx));
plot(I(t(idx(2:end))), frequencies, 'LineWidth', 2);
plot(I(t(idx(2:end))), pi./sqrt(I(t(idx(2:end)))), 'LineWidth', 2);
xlabel('$I$','Interpreter','latex', 'FontSize', labelfont)
ylabel('$T$','Interpreter','latex', 'FontSize', labelfont)
legend('$\hat{T}$','$\frac{\pi}{\sqrt{I}}$', 'Interpreter','latex', 'FontSize', labelfont)

print(fI, '../Figures/ThetaNeuronResponseToCurrentPeriod.png', '-dpng', '-r300')

%% Functions:
function I = excitabilitycurrent(t)
    I = max(t/10, power(t,2));
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

function I = fI_current(t)
    I = t;
end
