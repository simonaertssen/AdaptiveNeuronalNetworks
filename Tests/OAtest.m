close all; clc;
%% Setup
addpath('../Functions');
addpath('../Mean Field Reductions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 15;
labelfont = 13;
export = true;

%% Testing the OA approach:
tnow = 0; tend = 10;
h = 0.005;

pars.N = 5000;
pars.a_n = 0.666666666666666666667;
pars.eta0 = 0.5; pars.delta = 0.7; pars.K = 2;

seed = 2; rng(seed);
IC = - pi/2 * ones(pars.N, 1);

sfpars = make_scalefreeparameters(pars, 3);
sfpars = prepareOAparameters(sfpars);

% The OA mean field theory:
options = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);

[TOA, ZOA] = OA_simulatenetwork2(tnow, tend, -1i, sfpars);


%% Functions
function p = prepareOAparameters2(p)
    p.k = unique(p.degrees_in);
    p.l = numel(p.k);
    pkperm = p.k(randperm(p.l));
    p.OA = zeros(p.l, p.l);
    for i = 1:p.l
        ks = p.k(i)*ones(p.l,1);
        p.OA(i, :) = p.P(p.k).*assortativity(p.k, pkperm, ks, ks, p.N, p.meandegree, 0)/p.meandegree;
    end
end

function dzdt = MFROA2(t, z, p)
% Here we compute the differential equation for the mean field reduction 
% using the formulation for different types of networks, Ott-Antonsen 2017
    one = -1i.*(z-1).*(z-1)/2;
    two = (z+1).*(z+1);
    zc = conj(z);
    H = (1 + (z.*z + zc.*zc)/6 - 4.*real(z)/3);
    HOA = p.OA*H;
    dzdt = one + two.*(-p.delta + 1i*p.eta0 + 1i*p.K.*HOA)/2;
end

function [TOA, ZOA, b] = OA_simulatenetwork2(tnow, tend, IC, p, odeoptions)
    if nargin < 5
        odeoptions = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
    end

    if numel(IC) > 1
        OAIC = zeros(1,p.l);
        for i = 1:p.l
            OAIC(i) = sum(exp(1i*IC(p.degrees_in == p.k(i)))) / p.P(p.k(i));
        end
    elseif numel(IC) == 1
        OAIC = IC*ones(1,p.l);
    else
        error('IC might be wrong?')
    end
    
    [TOA, b] = ode45(@(t,x) MFROA2(t,x,p), [tnow, tend], gather(OAIC), odeoptions);
    ZOA = b*p.P(p.k)/p.N;
end

