clear all; close all; clc;
% In this script we will test a 2D version of the OA jacobian, to gain
% insight on the stability of fixpoints and to actually find them.

%% Setup:
addpath('../Functions');
addpath('../Mean Field Reductions/');

pars.N = 5000;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
pars.eta0 = 0.4; pars.delta = 0.7; pars.K = 2;
% pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;

seed = 1; rng(seed);
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
p = prepareOAparameters(make_randomparameters(pars, 0.33));

%%

figure; hold on
phasespaceplot();

Z0 = 0.3 - 1i*0.5;
zoa = map_Ztozoa(Z0, p)';
scatter(real(Z0), imag(Z0), 150, '+k')
scatter(real(zoa*p.P(p.k)/p.N), imag(zoa*p.P(p.k)/p.N), 150, 'xb');

eqpts = OA_fixedpointiteration(zoa', p);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 150, '+r');

% eqpts = NewtonRaphsonIteration(zoa', p);
% ZOA = eqpts'*p.P(p.k)/p.N;
% scatter(real(ZOA), imag(ZOA), 150, 'xb');

[TOA, b] = ode45(@(t,z) MFROA2Dtest(t,z,p), [0, 10], zoa);
Z = b*p.P(p.k)/p.N;
scatter(real(Z(1)), imag(Z(1)), 100, 'ob');
plot(real(Z), imag(Z));


%%
function J = MFROAJ(z,p)
    x = real(z); y = imag(z);
    dzdt = zeros(p.Mk, p.Mk);
    
end

function [x, xs] = NewtonRaphsonIteration(x0, p)

    f = @(z, p) MFROA2Dtest(0,z,p);
    
    function dfdz = df(z, p)
        zc = conj(z);
        H = (1 + (z.*z + zc.*zc)/6 - 4.*real(z)/3);
        I = -p.delta + 1i*p.eta0 + 1i*p.OA*H;
        dfdz = zeros(p.Mk, p.Mk);
        
        for r = 1:p.Mk
            zr = z(r);
            for c = 1:p.Mk
                dfdz(r,c) = 0.5*(z(r)+1)^2 * (1i*p.OA(r,c)*(z(c)-2)/3);
%                 Itmp = 1i*p.K/p.meandegree*(p.P(p.k(c))*assortativity(p.k(c), p.k_o(c), p.k(r), p.k_o(r), p.N, p.meandegree, 0)*(z(c)-2)/3);
%                 dfdz(r,c) = 0.5*(zr+1)^2 * Itmp;
                if r == c
                    dfdz(r,c) = dfdz(r,c) - 1i*(zr-1) + (zr+1)*I(r);
                end
            end
        end
    end
    
    x = x0;
    maxevals = 50;
    xs = x0'*p.P(p.k)/p.N;
    for evaltime = 1:maxevals
        x0 = x;
        x = x - df(x,p)\f(x,p);

        error = norm(x - x0);
        if error < 1.0e-10
            break
        end
        xs = [xs, x'*p.P(p.k)/p.N];
    end
    plot(real(xs), imag(xs), 'LineWidth', 2)
    disp(['Algorithm took ', num2str(evaltime), ' steps'])
    test = df(x0, p);
    test(1:5,1:5)
end