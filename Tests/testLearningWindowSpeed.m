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

%% Observe function shape
dt = linspace(-0.08, 0.1, tpts);

f_windows = figure('Renderer', 'painters', 'Position', [50 800 800 400]); box on; hold on;
ylim([-0.14, 0.12]);

plot(dt, Kempter1999Window(dt), 'LineWidth', 2, 'Color', '#77AC30')

yleft = Song2000Window(dt(dt<=0)); yright = Song2000Window(dt(dt>0));
plot(dt(dt<=0), yleft, 'LineWidth', 2, 'Color', '#0072BD'); 
plot(dt(dt>0), yright, 'LineWidth', 2, 'Color', '#0072BD', 'HandleVisibility', 'off')
line([0, 0],[min(yleft), max(yright)],'Color','k','LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off')

y = ChrolCannon2012Window(dt);
plot(dt, y, 'LineWidth', 2, 'Color', '#D95319')

plot(dt, Waddington2014Window(dt), 'LineWidth', 2, 'Color', '#77AC30')

legend("$$W(t)_S$$", "$$W(t)_C$$", "$$W(t)_W$$", 'Interpreter', 'latex', 'FontSize', labelfont)
xlabel("$t$ [s]", 'Interpreter', 'latex', 'FontSize', labelfont); 
ylabel("$W$", 'Interpreter', 'latex', 'FontSize', labelfont);
% exportpdf(f_windows, '../Figures/LearningWindows.pdf', export);
% close(f_windows)

%% Function properties
clc
% integral(@Kempter1999Window, -100, 100)
integral(@Song2000Window, -0.08, 0.1)
integral(@ChrolCannon2012Window, -0.08, 0.1)
integral(@Waddington2014Window, -0.08, 0.1)

%% Evaluate multiple time:
clc
tic;
timeme(N, tpts, @Waddington2014Window);
toc;

tic;
timeme(N, tpts, @ChrolCannon2012Window);
toc;

tic;
timeme(N, tpts, @Kempter1999Window);
toc;

% Turns out Waddington is quicker

%% Functions
function out = timeme(N, tpts, handle)
    dt = randn(N,1);
    for i = 1:tpts
        out = handle(dt);
    end
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
