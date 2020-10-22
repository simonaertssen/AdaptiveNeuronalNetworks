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

pars.N = 50;
pars.a_n = 0.666666666666666666667;
pars.eta0 = 0.5; pars.delta = 0.7; pars.K = 2;

seed = 2; rng(seed);
IC = - pi/2 * ones(pars.N, 1);

% Old version:
sfpars = make_scalefreeparameters(pars, 3);
sfpars = prepareOAparameters(sfpars);

% The OA mean field theory:
[TOA, ZOA] = OA_simulatenetwork(tnow, tend, -1i, sfpars);
plot(TOA, abs(ZOA), 'k')

% New version: a simulation per (k_in, k_out)
sfpars = make_scalefreeparameters(pars, 3);
sfpars = prepareOAparameters2(sfpars);

% The OA mean field theory:
[TOA, ZOA] = OA_simulatenetwork2(tnow, tend, -1i, sfpars);
plot(TOA, abs(ZOA), 'r')


%% Functions
pars.N = 100
sfpars = make_scalefreeparameters(pars, 3);
[TOA, ZOA] = OA_simulatenetwork2(tnow, tend, -i, sfpars);

function [TOA, ZOA, b] = OA_simulatenetwork2(tnow, tend, IC, p, odeoptions)
    if nargin < 5
        odeoptions = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
    end
        
    p.l_i = numel(unique(p.degrees_i));
    p.l_o = numel(unique(p.degrees_o));
    p.k = unique([p.degrees_i, p.degrees_o], 'rows');
    [p.l, ~] = size(p.k);
    p.l
    
    if numel(IC) > 1
        OAIC = zeros(1,p.l);
        for i = 1:p.l
            idx = (p.degrees_i == p.k(i, 1) & p.degrees_o == p.k(i, 2));
            size(idx)
            sum(exp(1i*IC(idx)))
            OAIC(i) = sum(exp(1i*IC(idx)) / (p.P(p.k(i,1)) * p.P(p.k(i,2))));
        end
    elseif numel(IC) == 1
        OAIC = IC*ones(1,p.l);
    else
        error('IC might be wrong?')
    end
    size(IC)
%     [TOA, b] = ode45(@(t,x) MFROA2(t,x,p), [tnow, tend], gather(OAIC), odeoptions);
%     ZOA = b*p.P(p.k)/p.N;
end

function dzdt = MFROA2(t, z, p)
% Here we compute the differential equation for the mean field reduction 
% using the formulation for different types of networks, Ott-Antonsen 2017
    one = -1i.*(z-1).*(z-1)/2;
    two = (z+1).*(z+1);
    zc = conj(z);
    H = (1 + (z.*z + zc.*zc)/6 - 4.*real(z)/3);
    
%     HOA = p.OA*H;
    HOA = zeros(p.l)
    for i = 1:p.l_i
        
    end
    dzdt = one + two.*(-p.delta + 1i*p.eta0 + 1i*p.K.*HOA)/2;
end


function p = prepareOAparameters2(p)
    p.pairs = unique([p.degrees_i, p.degrees_o], 'rows');
    p.l = numel(p.pairs); 
    p.k_i = unique(p.degrees_i);
    p.k_o = unique(p.degrees_o);
    
    products = zeros(size(p.pairs));
    for i = 1:p.l
        products(i,1) = p.P(p.k_i)*assortativity(p.k_i, pkperm, ks, ks, p.N, p.meandegree, 0)
        products(i,2) = p.P(p.k_i)*assortativity(p.k_i, pkperm, ks, ks, p.N, p.meandegree, 0)
    end
    p.products = products
    
    p.l_i = numel(p.k);
    pkperm = p.k(randperm(p.l_i));
    p.OA = zeros(p.l_i, p.l_i);
    for i = 1:p.l
        ks = p.k(i)*ones(p.l,1);
        p.OA(i, :) = p.P(p.k).*assortativity(p.k, pkperm, ks, ks, p.N, p.meandegree, 0)/p.meandegree;
    end
end



% Deprecated, only for inspiration
function scalefreepars = make_scalefreeparameters2(pars, shouldbesame, degree_i, degree_o, krange_i, krange_o)
    scalefreepars = pars;
    fsolveoptions = optimset('Display','off');
    % Make the necessery parameters for the scalefree networks
    if nargin < 1; error('Not enough input arguments'); end
    if nargin < 2; shouldbesame = false; end
    if nargin < 3; degree_i = 3; end
    if nargin < 4; degree_o = degree_i; end
    if nargin < 5; krange_i = [round(pars.N*3/20); round(pars.N*2/5)]; end
    if nargin < 6; krange_o = krange_i; end
   
    if degree_i < 2 || degree_o < 2
        error('Scale free networks do not exist for degress less than 2');
    end
    
    if degree_i == degree_o && krange_i(1) == krange_o(1) && krange_i(2) == krange_o(2)
        allequal = true;
    else
        allequal = false;
    end
    
    scalefreepars.krange_i = krange_i;
    scalefreepars.krange_o = krange_o;
    scalefreepars.degree_i = degree_i;
    scalefreepars.degree_o = degree_i;
    
    scalefreepars.P_i = @(x) scalefreepdf(x, pars.N, degree_i, krange_i(1), krange_i(2));
    scalefreepars.degrees_i = randsample(krange_i(1):krange_i(2), pars.N, true, scalefreepars.P_i(krange_i(1):krange_i(2)))';
    
    if allequal
        scalefreepars.P_o = scalefreepars.P_i;
        scalefreepars.degrees_o = scalefreepars.degrees_i(randperm(pars.N));
    else
        scalefreepars.P_o = @(x) scalefreepdf(x, pars.N, degree_o, krange_o(1), krange_o(2));
        scalefreepars.degrees_o = randsample(krange_o(1):krange_o(2), pars.N, true, scalefreepars.P_o(krange_o(1):krange_o(2)))';
    end
    
    % Assert maximum degree
    if max(scalefreepars.degrees_i) > pars.N - 1
        disp(['Setting higher in-degrees to ', num2str(pars.N-1)]);
        scalefreepars.degrees_i(scalefreepars.degrees_i > pars.N - 1) = pars.N - 1;
    end
    if max(scalefreepars.degrees_o) > pars.N - 1
        disp(['Setting higher in-degrees to ', num2str(pars.N-1)]);
        scalefreepars.degrees_o(scalefreepars.degrees_o > pars.N - 1) = pars.N - 1;
    end
    
    % Assert numer of links is equal
    if shouldbesame == true && sum(scalefreepars.degrees_i) ~= sum(scalefreepars.degrees_o)
        d_i_n = sum(scalefreepars.degrees_i);
        d_o_n = sum(scalefreepars.degrees_o);
        if d_i_n > d_o_n
            ntoclear = d_i_n - d_o_n
            if fix(ntoclear/pars.N) > 0
                scalefreepars.degrees_i = scalefreepars.degrees_i - fix(ntoclear/pars.N);
            end
            idx = randperm(pars.N,rem(ntoclear,pars.N));
            scalefreepars.degrees_i(idx) = scalefreepars.degrees_i(idx) - 1;
        else
            ntoclear = d_o_n - d_i_n
            if fix(ntoclear/pars.N) > 0
                scalefreepars.degrees_o = scalefreepars.degrees_o - fix(ntoclear/pars.N);
            end
            idx = randperm(pars.N,rem(ntoclear,pars.N));
            scalefreepars.degrees_o(idx) = scalefreepars.degrees_o(idx) - 1;
        end
    end
    
    scalefreepars.meandegree_i = fsolve(@(z) scalefreepars.P_i(z) - mean(scalefreepars.P_i(krange_i(1):krange_i(2))), krange_i(1), fsolveoptions);
    if allequal
        scalefreepars.meandegree_o = scalefreepars.meandegree_i;
    else
        scalefreepars.meandegree_o = fsolve(@(z) scalefreepars.P_o(z) - mean(scalefreepars.P_o(krange_o(1):krange_o(2))), krange_o(1), fsolveoptions);
    end
end

