close all; clear all; clc;
addpath('../Functions');

%% Test the different networks and their statistical properties
% We need to evaluate whether some of the random numbers we are pulling are
% actually following the right pdf.

pars.N = 10000;

%% A fixed degree network / diracnet: CORRECT
fixeddegreepars = make_fixeddegreeparameters(pars, round(pars.N/2));

if sum(fixeddegreepars.P(1:pars.N)) == pars.N; disp('Sum is correct'); end
if fixeddegreepars.meandegree == sum(fixeddegreepars.degrees_i)/pars.N; disp('Mean degree is correct'); end

limits = round(fixeddegreepars.meandegree + fixeddegreepars.meandegree*[-0.002, 0.002]);

figure; hold on; grid on;
xlim(limits);
x = limits(1):limits(2);
histogram(fixeddegreepars.degrees_i, 'Normalization', 'pdf');
plot(x, fixeddegreepars.P(x)/pars.N);
xline(fixeddegreepars.meandegree, 'k--')

%% A random network: CORRECT
randompars = make_randomparameters(pars, 0.02);

if round(sum(randompars.P(1:pars.N))) == pars.N; disp('Sum is correct'); end
if abs(randompars.meandegree - sum(randompars.degrees_i)/pars.N) < 1.0e-3*randompars.meandegree; disp('Mean degree is correct'); end

limits = round(randompars.meandegree + randompars.meandegree*[-1, 1]);

figure; hold on; grid on;
xlim(limits);
x = limits(1):limits(2);
histogram(randompars.degrees_i, 'Normalization', 'pdf');
plot(x, randompars.P(x)/pars.N);
xline(randompars.meandegree, 'k--')

%% A scale free network: CORRECT
kmin = 1000; kmax = 2000;
scalefreepars = make_scalefreeparameters(pars, 3, kmin, kmax);

if round(sum(scalefreepars.P(1:pars.N))) == pars.N; disp('Sum is correct'); end
if abs(scalefreepars.meandegree - sum(scalefreepars.degrees_i)/pars.N) < 1.0e-3*scalefreepars.meandegree; disp('Mean degree is correct'); end

limits = [kmin, kmax];

figure; hold on; grid on;
x = limits(1):limits(2);
histogram(scalefreepars.degrees_i, 'Normalization', 'pdf');
plot(x, scalefreepars.P(x)/pars.N);
xline(scalefreepars.meandegree, 'k--')

%% Test: a lognormal network - inspired by scale-free CORRECT
lognormpars = make_lognormparameters(pars, 3, 1, 500);

if round(sum(lognormpars.P(1:pars.N))) == pars.N; disp('Sum is correct'); end
if abs(lognormpars.meandegree - sum(lognormpars.degrees_i)/pars.N) < 1.0e-3*lognormpars.meandegree; disp('Mean degree is correct'); end

limits = [lognormpars.kmin-10, lognormpars.kmin+200];

figure; hold on; grid on;
x = limits(1):0.01:limits(2);
xlim([limits(1), limits(2)]);
histogram(lognormpars.degrees_i, 'Normalization', 'pdf');
plot(x, lognormpars.P(x)/pars.N);
xline(lognormpars.meandegree, 'k--')
