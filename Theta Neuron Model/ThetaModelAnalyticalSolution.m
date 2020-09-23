% We will test whether the analytical solution 
% theta = 2*np.arctan(np.sqrt(i)*np.tan(t*np.sqrt(i) - np.pi/2))
% is valid by recording the error between simullations and theory.

clear all; close all; clc;

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

tnow = 0; tend = 50;
F = @thetaneuron; h = 0.001; IC = -pi;

fanalytic = figure('Renderer', 'painters', 'Position', [50 800 1000 400]);
titlefont = 15;
labelfont = 15;
m = 2; n = 1;


%% A simple constant current

I = @simplecurrent;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, IC, h, I);
drawthetas = spikesNaN(thetas);
realthetas = solution(t,I);

imrow(1) = subplot(m,n,1); hold on; box on;

yyaxis left
ylim([-pi - 1.0, pi + 0.3]);
plot(t, thetas, ':k', 'LineWidth', 1, 'HandleVisibility','off');
plot(t, drawthetas, '-', 'LineWidth', 2, 'color', '#0072BD');

plot(t, realthetas, ':k', 'LineWidth', 1, 'HandleVisibility','off');
plot(t, spikesNaN(realthetas), '-', 'LineWidth', 2, 'color', '#77AC30');

ylabel('$\theta$','Interpreter','latex', 'FontSize', labelfont);
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)

yyaxis right
maxy = max(I(t));
plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
ylabel('$I$','Interpreter','latex', 'FontSize', labelfont);
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';
set(gca,'YTick', 0:maxy:maxy, 'YLim', [0, maxy*10]);

legend('Simulation', 'Analytical');


%% A more difficult sinusoidal current

I = @sinecurrent;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, IC, h, I);
drawthetas = spikesNaN(thetas);
realthetas = solution(t,I);

imrow(2) = subplot(m,n,2); hold on; box on;

yyaxis left
ylim([-pi - 1.0, pi + 0.3]);
plot(t, thetas, ':k', 'LineWidth', 1, 'HandleVisibility','off');
plot(t, drawthetas, '-', 'LineWidth', 2, 'color', '#0072BD');

plot(t, realthetas, ':k', 'LineWidth', 1, 'HandleVisibility','off');
plot(t, spikesNaN(realthetas), '-', 'LineWidth', 2, 'color', '#77AC30');

ylabel('$\theta$','Interpreter','latex', 'FontSize', labelfont);
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)

yyaxis right
maxy = max(I(t));
plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
ylabel('$I$','Interpreter','latex', 'FontSize', labelfont);
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';
set(gca,'YTick', 0:maxy:maxy, 'YLim', [0, maxy*10]);

legend('Simulation', 'Analytical');


%% Save the figure:
print(fanalytic, '../Figures/ThetaNeuronAnalyticalSolution.png', '-dpng', '-r300')


%% Functions
function realthetas = solution(t,I)
    s = sqrt(I(t));
    realthetas = 2*atan(-s.*cot(t.*s));
end


function I = simplecurrent(t)
    I = 2*ones(size(t));
end

function I = sinecurrent(t)
    I = 1.1 + sin(t);
end