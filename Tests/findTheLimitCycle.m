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
oscillation1 = sin(2*pi*t/11);
oscillation2 = 0.1*sin(2*pi*t/5);
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
for i = 1:2
    if sign(dfdetr(t_fmean(end))) == sign(dfdetr(t_fmean(end-i)))
        tstart = t_fmean(end-i)
    end  
end

% Test:
ts = findlimitcycle(func)
hold on;
plot(t, func)
plot(t(ts(1):ts(2)), func(ts(1):ts(2)))
hold off

%% To the frequency domain:
Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:n-1)*T;      
Y = fft(func);

P2 = abs(Y/n);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(n/2))/n;

plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Get the highest waves out:
clc
[peakvals, indices] = findpeaks(P1,f, 'NPeaks', 5)
%[peakvals, indices] = maxk(P1,5)
1./f(indices)


function ts = findlimitcycle(vector)
    vector = detrend(vector, 1);
    ts = [nan, nan]
    t_fmean = find(diff(sign(vector))); 
    dvector = diff(vector);

    tstart = nan;
    for i = 1:2
        if sign(dvector(t_fmean(end))) == sign(dvector(t_fmean(end-i)))
            tstart = t_fmean(end-i);
        end  
    end
    ts = [tstart, t_fmean(end)];
end