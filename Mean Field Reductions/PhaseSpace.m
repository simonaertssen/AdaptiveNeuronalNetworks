clear all; close all; clc;
% Here we will make phase space plots for the report, illustrating all the
% different global states

%%
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')
rect = [-50 50 800 800];
cm = [1,0,0; 0, 0.7410, 0.4470; 0, 0.4470, 0.6410];

titlefont = 15;
labelfont = 13;
export = true;

%% Grids:
% Axis square:
unitsquare = [-1, -1, 1, 1; 1, -1, -1, 1];
th = 0:pi/50:2*pi;
unitcircle = [cos(th); sin(th)];
drawcircle = round(unitcircle);

% Main grid:
l = 1; stp = 2*l/40; interval = -l:stp:l;
[X,Y] = meshgrid(interval,interval);

[in, ~] = inpolygon(X, Y, cos(th), sin(th));
[sz, ~] = size(in); m = ceil(sz/2);
in(1, :) = 0; in(end, :) = 0; in(:, 1) = 0; in(:, end) = 0;
in(1, m) = 1; in(end, m) = 1; in(m, 1) = 1; in(m, end) = 1;

% Draw grid
l = 1; stp = 2*l/6; interval = -l:stp:l;
[X1,Y1] = meshgrid(interval,interval);

% Startpoints
[in, on] = inpolygon(X1, Y1, [0.9*cos(th) 0], [0 0.9*sin(th)]);
[sz, ~] = size(in); m = ceil(sz/2);
in(1, :) = 0; in(end, :) = 0; in(:, 1) = 0; in(:, end) = 0;
in(1, m) = 1; in(end, m) = 1; in(m, 1) = 1; in(m, end) = 1;
X1start = reshape(X1(in), 5, 5); Y1start = reshape(Y1(in), 5, 5);

%% Theta neurons parameters:
pars.N = 1000;
pars.a_n = 0.666666666666666666667;
seed = 2; rng(seed);

%% 1. MFR phase space: PSR
pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

close all;
f_MFRPSR = figure('Renderer', 'painters', 'Position', rect); hold on; box on;

% Equilibrium point:
fsolveoptions = optimset('Display','off');
[eqpt,~,~,~,JACOB] = fsolve(@(x) MFR2D(0,x,pars), [-0.5 -0.5], fsolveoptions);

% Arrows:
complex = MFR(0, X + 1i*Y, pars);
quiver(X, Y, real(complex), imag(complex), 'color', cm(2,:))

% Lines:
complex = MFR(0, X + 1i*Y, pars);
startx = 0.8*cos( -pi/5:pi/5:pi); starty = 0.8*sin(-pi/5:pi/5:pi);
startx(2) = []; starty(2) = [];
scatter(startx, starty, 50, cm(3,:), 'filled', 'o', 'LineWidth',2);
s = streamline(X, Y, real(complex), imag(complex), startx, starty, [0.05,1000]);
set(s, 'color', cm(3,:), 'LineWidth', 2);
for i = 1:6
    dist = (s(i).XData - eqpt(1)).^2 + (s(i).YData - eqpt(2)).^2 - 0.02;
    dist(end) = -1;
    idx = find(dist < 0, 1, 'first');
    s(i).XData = s(i).XData(1:idx);
    s(i).YData = s(i).YData(1:idx);
    plot_arrow(s(i).XData(end-8), s(i).YData(end-8), s(i).XData(end), s(i).YData(end),'linewidth', 2, ...
    'color', cm(3,:),'facecolor', cm(3,:),'edgecolor', cm(3,:), 'headwidth',0.5,'headheight',1.5);
end

phasespaceplot();

% Plot the equilibrium points and Jacobian:
[V,D] = eig(JACOB);
index = 1;
for i = 1:2
    lambda = D(i,i);
    bef = 0.04;
    endl = 0.075;
    draw = 0.75*endl*lambda;
    if lambda < 0
        plot_arrow(eqpt(1) + endl*V(1,index), eqpt(2) + endl*V(2,index), eqpt(1) + bef*V(1,index), eqpt(2) + bef*V(2,index), 'headwidth',0.6,'headheight',1.5); 
        plot_arrow(eqpt(1) - endl*V(1,index), eqpt(2) - endl*V(2,index), eqpt(1) - bef*V(1,index), eqpt(2) - bef*V(2,index), 'headwidth',0.6,'headheight',1.5); 
        plot([eqpt(1) - draw*V(1,index), eqpt(1) + draw*V(1,index)], [eqpt(2) - draw*V(2,index), eqpt(2) + draw*V(2,index)], 'k', 'LineWidth', 2);
    end
    index = index + 1;
end
scatter(eqpt(1), eqpt(2), 150, 'or', 'filled')

% End figure:
set(gcf,'color','w'); set(gca,'FontSize',14); xlim([-1,1]); ylim([-1,1]);
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', 20)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', 20)
print(f_MFRPSR, '../Figures/MFR_PSR.png', '-dpng', '-r300')
close(f_MFRPSR)

%% 2. MFR phase space: PSS
pars.eta0 = 0.5; pars.delta = 0.7; pars.K = 2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

close all
f_MFRPSS = figure('Renderer', 'painters', 'Position', rect); hold on; box on;

% Equilibrium point:
fsolveoptions = optimset('Display','off');
[eqpt,~,~,~,JACOB] = fsolve(@(x) MFR2D(0,x,pars), [-0.5 -0.5], fsolveoptions);

% Arrows:
complex = MFR(0, X + 1i*Y, pars);
quiver(X, Y, real(complex), imag(complex), 'color', cm(2,:))

% Lines:
complex = MFR(0, X + 1i*Y, pars);
startx = 1; starty = 0;
scatter(startx, starty, 50, cm(3,:), 'filled', 'o', 'LineWidth',2);% 'color', cm(3,:));
s = streamline(X, Y, real(complex), imag(complex), startx, starty, [0.05,1550]);
set(s, 'color', cm(3,:), 'LineWidth', 2);
plot_arrow(s.XData(end-8), s.YData(end-8), s.XData(end), s.YData(end),'linewidth', 2, ...
    'color', cm(3,:),'facecolor', cm(3,:),'edgecolor', cm(3,:), 'headwidth', 0.5,'headheight',3);

phasespaceplot();

% Plot the equilibrium points and Jacobian:
[V,D] = eig(JACOB);
scatter(eqpt(1), eqpt(2), 150, 'or', 'filled')

% End figure:
hold off; set(gcf,'color','w'); set(gca,'FontSize',14); xlim([-1,1]); ylim([-1,1]); axis square;
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', 20)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', 20)
print(f_MFRPSS, '../Figures/MFR_PSS.png', '-dpng', '-r300')
close(f_MFRPSS)

%% 3. MFR phase space: CPW
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

% Start figure:
close all
f1 = figure('Renderer', 'painters', 'Position', rect);
% title('PSR state','Interpreter','none', 'FontSize', 20, 'FontWeight', 'normal')

hold on; box on;

% Equilibrium point:
fsolveoptions = optimset('Display','off');
[eqpt,FVAL,EXITFLAG,OUTPUT,JACOB] = fsolve(@(x) MFR2D(0,x,pars), [-0.8 -0.6], fsolveoptions);

% Arrows:
complex = MFR(0, X + 1i*Y, pars);
quiver(X, Y, real(complex), imag(complex), 'color', cm(2,:))

% Lines:
complex = MFR(0, X + 1i*Y, pars);
startx = [0.34, 0, 0, -0.6, -0.8]; starty = [-0.2, -1, -0.8, 0.4, 0.2];
scatter(startx, starty, 50, cm(3,:), 'filled', 'o', 'LineWidth',2);% 'color', cm(3,:));
s = streamline(X, Y, real(complex), imag(complex), startx, starty, [0.01,9000]);
set(s, 'color', cm(3,:), 'LineWidth', 2);

% Streamlines yield inexact results?
for i = 1:5
%     [TT, Z] = ode45(@(t,x) MFR(t,x,pars), [0, 5], startx + starty*1i, options);
%     disp("Fully connected network done")
%     plot(real(Z), imag(Z), 'k');
    
    dist = (s(i).XData - eqpt(1)).^2 + (s(i).YData - eqpt(2)).^2 - 0.02;
    dist(end) = -1;
    idx = find(dist < 0, 1, 'first');
    s(i).XData = s(i).XData(1:idx);
    s(i).YData = s(i).YData(1:idx);
    plot_arrow(s(i).XData(end-40), s(i).YData(end-40), s(i).XData(end), s(i).YData(end),'linewidth', 2, ...
    'color', cm(3,:),'facecolor', cm(3,:),'edgecolor', cm(3,:), 'headwidth',0.5,'headheight',1.5);     
end

% Fill the outside
fill([drawcircle(1,:) flip(unitcircle(1,:))], [drawcircle(2,:) flip(unitcircle(2,:))], 'w')

% Imaginary unit circle:
th = 0:pi/50:2*pi;
circ = plot(cos(th), sin(th), '-k', 'LineWidth', 2);

% Plot the equilibrium points and Jacobian:
ICs = [-0.5 -0.5; -0.8 -0.6];
for eqpointindex = 1:2
    [eqpt,FVAL,EXITFLAG,OUTPUT,JACOB] = fsolve(@(x) MFR2D(0,x,pars), ICs(eqpointindex,:), fsolveoptions);
    [V,D] = eig(JACOB);
    index = 1;
    for i = 1:2
        lambda = D(i,i);
        bef = 0.04;
        endl = 0.075;
        draw = 0.5*endl*lambda;
        if lambda < 0
            plot_arrow(eqpt(1) + endl*V(1,index), eqpt(2) + endl*V(2,index), eqpt(1) + bef*V(1,index), eqpt(2) + bef*V(2,index), 'headwidth',0.6,'headheight',1.5);
            plot_arrow(eqpt(1) - endl*V(1,index), eqpt(2) - endl*V(2,index), eqpt(1) - bef*V(1,index), eqpt(2) - bef*V(2,index), 'headwidth',0.6,'headheight',1.5);
            plot([eqpt(1) - draw*V(1,index), eqpt(1) + draw*V(1,index)], [eqpt(2) - draw*V(2,index), eqpt(2) + draw*V(2,index)], 'k', 'LineWidth', 2);
        else
            plot_arrow(eqpt(1) + endl*V(1,index), eqpt(2) + endl*V(2,index),eqpt(1) + draw*V(1,index), eqpt(2) + draw*V(2,index),'headwidth',0.6,'headheight',1.5);
            plot_arrow(eqpt(1) - endl*V(1,index), eqpt(2) - endl*V(2,index),eqpt(1) - draw*V(1,index), eqpt(2) - draw*V(2,index),'headwidth',0.6,'headheight',1.5);
            plot([eqpt(1) - draw*V(1,index), eqpt(1) + draw*V(1,index)], [eqpt(2) - draw*V(2,index), eqpt(2) + draw*V(2,index)], 'k', 'LineWidth', 2);
        end
        index = index + 1;
    end
    scatter(eqpt(1), eqpt(2), 150, 'or', 'filled')
end

% Draw central focus and draw round vector
[eqpt,FVAL,EXITFLAG,OUTPUT,JACOB] = fsolve(@(x) MFR2D(0,x,pars), [0,0], fsolveoptions);

t = linspace(0, 2*pi + 0.3, 100);
radius = t/max(t)/10;
curve = [eqpt(1) + radius.*cos(t); eqpt(2) + radius.*sin(t)];
plot(curve(1,:), curve(2,:), 'k', 'LineWidth', 2)

scatter(eqpt(1), eqpt(2), 150, 'or', 'filled')
plot_arrow(curve(1,end-6), curve(2,end-6), curve(1,end), curve(2,end), 'headwidth',0.6,'headheight',1.5);


% End figure:
hold off; set(gcf,'color','w'); set(gca,'FontSize',14); xlim([-1,1]); ylim([-1,1]); axis square;
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', 20)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', 20)
% text(-0.8, -0.8, ['$\eta = $' num2str(pars.eta0),], 'Interpreter','latex', 'HorizontalAlignment', 'center', 'FontSize',20);
% text(-0.8, -0.85, ['$\Delta = $' num2str(pars.delta),], 'Interpreter','latex', 'HorizontalAlignment', 'center', 'FontSize',20);
% text(-0.82, -0.9, ['$K = $' num2str(pars.K),], 'Interpreter','latex', 'HorizontalAlignment', 'center', 'FontSize',20);

% print(f1, 'figures/CPW.png', '-dpng', '-r400')

