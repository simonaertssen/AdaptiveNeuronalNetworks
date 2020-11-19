clear all; close all; clc;
% Here we will try to develop a method to accurately and quickly detect and
% find a limit cycle in a given vector, assuming it is real. We want to obtain 
% the last full revolution of the cycle.

%% Setup: making the corrupted signal
% The spine of this function is going to be a harmonic wave mad up of three
% distict oscillating frequencies, plus a bias that disappears over time.
% This emulates the behaviour of a real cycle pretty well.

n = 1000;
t = linspace(1,200,n);

inverse = 1./(t.^(0.2));
oscillation1 = sin(2*pi*t/7);
oscillation2 = 0.5*sin(2*pi*t/14);
oscillation3 = 0.01*sin(2*pi*t/111);

func = inverse + oscillation1 + oscillation2;
fdetr = detrend(func, 1);

cla
hold on;
plot(t, func)
yline(mean(func))
plot(t, fdetr)
yline(mean(fdetr))
hold off;

%% Take the last full revolution of the function:
% Start from the last time that the function went through the mean:
t_fmean = find(diff(sign(fdetr))); % last change = t_fmean(end)
dfdetr = diff(fdetr);

tend = t_fmean(end);
tstart = nan;
for i = 3:4
    if sign(dfdetr(t_fmean(end))) == sign(dfdetr(t_fmean(end-i)))
        tstart = t_fmean(end-i);
    end  
end

% Test:
ts = findlimitcycle(func)
tstart = ts(1); tend = ts(2);

cla
hold on;
plot(t, func)
plot(t(tstart:tend), func(tstart:tend))
hold off

%% To the frequency domain:
filtered = fitsine(func, t, 0.01);

cla
hold on;
plot(t, func)
plot(t, filtered)
hold off


function ts = findlimitcycle(vector)
    vector = detrend(vector, 1);
    ts = [nan, nan];
    t_fmean = find(diff(sign(vector))); 
    dvector = diff(vector);

    tstart = nan;
    for i = 3:4
        if sign(dvector(t_fmean(end))) == sign(dvector(t_fmean(end-i)))
            tstart = t_fmean(end-i);
        end  
    end
    ts = [tstart, t_fmean(end)];
end


function [result, peaks, troughs] = fitsine(y, t, eps)
% Fast fourier-transform
f = fft(y);
l = length(y);
p2 = abs(f/l);
p1 = p2(1:ceil(l/2+1));
p1(2:end-1) = 2*p1(2:end-1);
freq = (1/mean(diff(t)))*(0:ceil(l/2))/l;

% Find maximum amplitude and frequency
maxPeak1 = p1 == max(p1(2:end)); % disregard 0 frequency!
maxAmplitude1 = p1(maxPeak1);     % find maximum amplitude
maxFrequency1 = freq(maxPeak1);   % find maximum frequency

p1(maxPeak1) = [];
maxPeak2 = p1 == max(p1(2:end)); % disregard 0 frequency!
maxAmplitude2 = p1(maxPeak2);     % find maximum amplitude
maxFrequency2 = freq(maxPeak2);   % find maximum frequency

% Initialize guesses
p = [];
% Sine 1:
p(1) = mean(y);         % vertical shift
p(2) = maxAmplitude1;    % amplitude estimate
p(3) = maxFrequency1;    % phase estimate
p(4) = 0;               % phase shift (no guess)
p(5) = 0;               % trend (no guess)
% Sine 2:
p(6) = maxAmplitude2;    % amplitude estimate
p(7) = maxFrequency2;    % phase estimate
p(8) = pi/2;            % phase shift (no guess)

% Create model
f = @(p) p(1) + p(2)*sin(p(3)*2*pi*t+p(4)) + p(6)*sin(p(7)*2*pi*t+p(8)); % + p(5)*t;
ferror = @(p) sum((f(p) - y).^2);
% Nonlinear least squares
% If you have the Optimization toolbox, use [lsqcurvefit] instead!
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',1e-25);
[param,fval,exitflag,output] = fminsearch(ferror,p,options);

% Calculate result
result = f(p);

% Find peaks
peaks = abs(sin(param(3)*2*pi*t+param(4)) - 1) < eps;

% Find troughs
troughs = abs(sin(param(3)*2*pi*t+param(4)) + 1) < eps;

end
