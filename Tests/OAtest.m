close all; clc;
%% Setup
addpath('../Functions');
addpath('../Mean Field Reductions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 15;
labelfont = 13;
export = true;

%% Theta model parameters
tnow = 0; tend = 10;
h = 0.005;

pars.N = 100;
pars.a_n = 0.666666666666666666667;
pars.eta0 = 0.5; pars.delta = 0.7; pars.K = 2;

seed = 2; rng(seed);
IC = - pi/2 * ones(pars.N, 1);

%% Testing the OA approach:
figure; hold on

tic;
% Old version:
sfpars = prepareOAparameters(make_scalefreeparameters(pars, 3));
[TOA, ZOA] = OA_simulatenetwork(tnow, tend, -1i, sfpars);
plot(TOA, abs(ZOA), 'k')
toc;

%% 
N = 1000; l = 100;
k = randi([sfpars.kmin, sfpars.kmax], l, 2);
k = k ./ sum(k);
ks = k(2, 1)*ones(l,2);
test1 = assortativity(k(:,2), k(:,1), ks(:,2), ks(:,1), N, round(0.3*N), 0)'
test2 = assortativity2(k, ks, N, round(0.3*N), 0)'


%%
rng(2)
tic;
% New version: a simulation per (k_in, k_out)
sfpars = prepareOAparameters2(make_scalefreeparameters(pars, 3));
[TOA, ZOA] = OA_simulatenetwork2(tnow, tend, -1i, sfpars);
plot(TOA, abs(ZOA), 'r')
toc;

function p = prepareOAparameters2(p)
    [d_i, d_o] = meshgrid(unique(p.degrees_i), unique(p.degrees_o));
    p.k = [reshape(d_i, numel(d_i), 1), reshape(d_o, numel(d_o), 1)];
    p.l = numel(p.k)/2;
%     p.P = @(x) p.P(x)*sum(p.P(p.k(:,1)))/p.N;
    
%     p.k = unique(p.degrees_i);
%     p.l = numel(p.k);
%     pkperm = p.k(randperm(p.l));
    p.OA = zeros(p.l, p.l);
    for i = 1:p.l
        a = p.P(p.k(:,1)).*p.P(p.k(:,2)).*assortativity2(p.k, p.k(i,:).*ones(p.l,2), p.N, p.meandegree, 0)/p.meandegree;
        p.OA(i,:) = a;
%         ks = p.k(i,1)*ones(p.l,1);
%         p.OA(i, :) = p.P(p.k(:,1)).*assortativity(p.k(:,1), p.k(:,1), ks, ks, p.N, p.meandegree, 0)/p.meandegree;
    end
end


function [TOA, ZOA, b] = OA_simulatenetwork2(tnow, tend, IC, p, odeoptions)
    if nargin < 5
        odeoptions = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
    end
        
    Ps = p.P(p.k(:,1)) .* p.P(p.k(:,2));
    if numel(IC) > 1
        OAIC = zeros(p.l,1);
        for i = 1:p.l
            idx = (p.degrees_i == p.k(i, 1) & p.degrees_o == p.k(i, 2));
            OAIC(i) = sum(exp(1i*IC(idx)) / Ps(i));
        end
    elseif numel(IC) == 1
        OAIC = IC*ones(p.l,1);
    else
        error('IC might be wrong?')
    end
    [TOA, b] = ode45(@(t,x) MFROA2(t,x,p), [tnow, tend], gather(OAIC), odeoptions);
    ZOA = b*Ps/(p.N*p.N);
end

function a = assortativity2(k_accent, k, N, k_mean, c)
if c == 0
    a = max(0, min(1, (k_accent(:,2).*k(:,1)/(N*k_mean))));
    a = k_accent(:,2).*k(:,1)/(N*k_mean);
% else
%     a = max(0, min(1, (k_accent(2).*k(1)/(N*k_mean)) .* (1 + c*((k_accent_in - k_mean)./k_accent_out).*((k_out - k_mean)./k_in))));
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
    dzdt = one + 0.5*two.*(-p.delta + 1i*p.eta0 + 1i*p.K.*HOA);
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

