clear all; close all; clc
% Test the evaluation speed of some of the learning windows:
%% Setup
addpath('../Functions');

N = 1000;
tpts = 30001;

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 15;
labelfont = 13;
export = true;

%% Observe function shapes: biphasic windows
dt = linspace(-0.08, 0.1, tpts);

f_windows = figure('Renderer', 'painters', 'Position', [50 800 600 200]); box on; hold on;
ylim([-0.14, 0.12]);

plot(dt, Kempter1999Window(dt), 'LineWidth', 2, 'Color', '#0072BD')

yleft = Song2000Window(dt(dt<=0)); yright = Song2000Window(dt(dt>0));
plot(dt(dt<=0), yleft, 'LineWidth', 2, 'Color', '#D95319'); 
plot(dt(dt>0), yright, 'LineWidth', 2, 'Color', '#D95319', 'HandleVisibility', 'off')
line([0, 0],[min(yleft), max(yright)],'Color','k','LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off')
scatter(0, Song2000Window(0), 30, [0.8500 0.3250 0.0980], 'filled', 'o')

legend("$$W(t)_K$$", "$$W(t)_S$$", 'Interpreter', 'latex', 'FontSize', labelfont, 'Orientation','horizontal')
xlabel("$t$ [s]", 'Interpreter', 'latex', 'FontSize', labelfont); 
ylabel("$W$", 'Interpreter', 'latex', 'FontSize', labelfont);

exportpdf(f_windows, '../Figures/LearningWindowsBiphasic.pdf', export);
close(f_windows)

%% Observe function shapes: triphasic windows
dt = linspace(-0.08, 0.1, tpts);

f_windows = figure('Renderer', 'painters', 'Position', [50 800 600 200]); box on; hold on;
ylim([-0.14, 0.12]);

y = ChrolCannon2012Window(dt);
plot(dt, y, 'LineWidth', 2, 'Color', '#77AC30')

plot(dt, Waddington2014Window(dt), 'LineWidth', 2, 'Color', '#A2142F')

legend("$$W(t)_C$$", "$$W(t)_W$$", 'Interpreter', 'latex', 'FontSize', labelfont, 'Orientation','horizontal')
xlabel("$t$ [s]", 'Interpreter', 'latex', 'FontSize', labelfont); 
ylabel("$W$", 'Interpreter', 'latex', 'FontSize', labelfont);
exportpdf(f_windows, '../Figures/LearningWindowsTriphasic.pdf', export);
close(f_windows)

%%
dt = linspace(-0.08, 0.1, tpts);
plot(dt, Song2017Window(dt))

%% Function properties
clc
integral(@Kempter1999Window, -0.08, 0.1)
integral(@Song2000Window, -0.08, 0.1)
integral(@Song2017Window, -0.08, 0.1)
integral(@ChrolCannon2012Window, -0.08, 0.1)
integral(@Waddington2014Window, -0.08, 0.1)

%% Evaluate multiple time:
clc

tic;
timeme(N, 1000, @Kempter1999Window);
toc;

tic;
timeme(N, 1000, @Song2000Window);
toc;

tic;
timeme(N, 1000, @ChrolCannon2012Window);
toc;

tic;
timeme(N, 1000, @Waddington2014Window);
toc;

% Turns out Waddington is quicker

%% Intrinsic plasticity:
% Use the theory treated in Song2017
dt = linspace(-0.01, 0.25, tpts);

plot(dt, Song2017ip(dt))

%% Functions
function out = timeme(N, tpts, handle)
    dt = randn(N,1);
    for i = 1:tpts
        out = handle(dt);
    end
end

function dW = Kempter1999Window(dt)
    t_syn = 5;
    t_pos = 1;
    t_neg = 20;
    A_p = 1;
    A_n = -1;
    learning_rate = 5.0e-2;
    eps = 1.0e-9;
    dt = -dt * 1.0e3; % Convert to seconds
    
    dW = zeros(size(dt));
    t_neg_idx = dt <= 0;
    t_pos_idx = 0 < dt;
    dW(t_neg_idx) = exp(dt(t_neg_idx)/t_syn + eps).*(A_p*(1-dt(t_neg_idx)/t_pos) + A_n*(1-dt(t_neg_idx)/t_neg));
    dW(t_pos_idx) = A_p*exp(-dt(t_pos_idx)/t_pos + eps) + A_n*exp(-dt(t_pos_idx)/t_neg + eps);
    dW = learning_rate * dW;
end


function dW = Song2000Window(dt)
    t_pos = 20;
    t_neg = 20;
    A_pos = 0.1;
    A_neg = -0.12;

    dt = dt * 1.0e3; % Convert to seconds
    dW = zeros(size(dt));
    t_zero = dt == 0;
    t_neg_idx = dt <= 0;
    t_pos_idx = 0 < dt;
    dW(t_zero) = 0;
    dW(t_pos_idx) = A_pos*exp(-dt(t_pos_idx)/t_pos);
    dW(t_neg_idx) = A_neg*exp(dt(t_neg_idx)/t_neg);
end


function dW = Song2017Window(dt)
    t_pos = 20;
    t_neg = 20;
    A_pos = 0.005;
    A_neg = -0.00525;

    dt = dt * 1.0e3; % Convert to seconds
    dW = zeros(size(dt));
    t_pos_idx = dt > 0;
    t_neg_idx = dt <= 0;
    dW(t_pos_idx) = A_pos*exp(-dt(t_pos_idx)/t_pos);
    dW(t_neg_idx) = A_neg*exp(dt(t_neg_idx)/t_neg);
end


function dW = Waddington2014Window(dt)
%     lr = 0.005; alpha = 8.0e-6;
    lr = 0.1; alpha = 4.0;
    dt = dt * 1.0e3; % Convert to seconds
    dW = lr*(1 - ((dt-alpha).^2)./alpha^2).*exp(-abs(dt-alpha)./alpha);
end


function dW = ChrolCannon2012Window(dt)
%     Ap = 2.3e-4; 
%     Am = 1.5e-4;
%     dW = Ap*exp(-((dt + 15e-7).^2/2e-10)) - Am*exp(-((dt + 20e-7).^2/2e-9));

    Ap = 0.23; 
    Am = 0.15;
    t_pos = 200;
    t_neg = 2000;
    dt = dt * 1.0e3; % Convert to seconds
    dW = Ap*exp(-((dt - 15).^2/t_pos)) - Am*exp(-((dt - 20).^2/t_neg));
end


function phi = Song2017ip(dt)
    T_max = 90; %ms
    T_min = 110; %ms
    lr = 0.012;

    dt = dt * 1.0e3; % Convert to seconds
    phi = zeros(size(dt));
    t_neg_idx = dt < T_min;    
    t_pos_idx = dt > T_max;
    
    phi((T_min <= dt) & (dt <= T_max)) = 0;
    phi(t_neg_idx) = -lr*exp((T_min - dt(t_neg_idx))/T_min);
    phi(t_pos_idx) = lr*exp((dt(t_pos_idx) - T_max)/T_max);
end
