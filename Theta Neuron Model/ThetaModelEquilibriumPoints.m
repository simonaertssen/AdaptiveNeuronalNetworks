%% Equilibrium points of the theta model

clear all; close all; clc;

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

tnow = 0; tend = 4;
F = @thetaneuron; h = 0.001;

feqpts = figure('Renderer', 'painters', 'Position', [50 800 1000 200]); 
titlefont = 15;
labelfont = 20;
axesfont = 15;
m = 1; n = 3;

%% Steady current: 0.5

currents = [-0.5, -1, -1.2];
for i = 1:3
    current = currents(i);
    I = @(t) steadycurrent(t, current);
    
    [t, thetas] = DOPRI_singleneuron(F, tnow, tend, [-pi; - pi/2; 0; pi/2; pi], h, I);
    drawthetas = spikesNaN(thetas);

    imrow(i) = subplot(m,n,i); hold on; box on;

%     title(["I = " num2str(-current)], 'FontSize', titlefont, 'FontName', 'SansSerif');

    yyaxis left
    ylim([-pi, pi]);
    plot(t, drawthetas, '-', 'LineWidth', 2, 'color', '#0072BD');
    eqpt = -acos((I(t) + 1)./(1 - I(t)));
    plot(t, eqpt, '--k', 'LineWidth', 1);
    text(tend - 0.1, eqpt(1) + 0.4, num2str(eqpt(1)), 'HorizontalAlignment', 'right');
    
    plot(t, -eqpt, '--k', 'LineWidth', 1);
    text(tend - 0.1, -eqpt(1) + 0.4, num2str(-eqpt(1)), 'HorizontalAlignment', 'right');
    
    xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
    set(gca,'YTick',-pi:pi/2:pi, 'YTickLabel',{'-$$\pi$$','$$\frac{\pi}{2}$$','0','$$\frac{\pi}{2}$$','$$\pi$$'}, 'TickLabelInterpreter', 'latex', 'YLim', [-pi-0.2, pi + 0.2], 'FontSize', axesfont)
    
    if i == 1
        ylabel('$\theta$','Interpreter','latex', 'FontSize', labelfont);
    end

  
    yyaxis right
    maxy = max(I(t));
    plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
    ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';
    set(gca,'YTick', -1:0.5:0, 'YLim', [-1.5, 4]);
    
    if i == 3
        ylabel('$I$','Interpreter','latex', 'FontSize', labelfont);
    end
    
%     removewhitspace();

end

%% Save:
% print(feqpts, '../Figures/ThetaModelEquilibriumPoints.png', '-dpng', '-r300')
exportgraphics(feqpts,'../Figures/ThetaModelEquilibriumPoints.pdf')


%% Functions:
function I = steadycurrent(t, current)
    I = current*ones(size(t));
end
