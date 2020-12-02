clear all; close all; clc;
% In this script we will test a 2D version of the OA jacobian, to gain
% insight on the stability of fixpoints and to actually find them.

%% Setup:
addpath('../Functions');
addpath('../Mean Field Reductions/');

pars.N = 5000;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9; % CPW (limit cycle)
% pars.eta0 = 0.4; pars.delta = 0.7; pars.K = 2;    % PSS (stable spiral)
% pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;  % PRS (stable focus)

seed = 1; rng(seed);
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
p = prepareOAparameters(make_randomparameters(pars, 0.4));

opts = odeset('RelTol', 1.0e-9,'AbsTol', 1.0e-9);

%%
clc

p

figure; hold on
phasespaceplot();

Z0 = -0.3 + 1i*0.9;
% zoa = map_Ztozoa(Z0, p)';
zoa = Z0*ones(1,p.Mk);
scatter(real(Z0), imag(Z0), 150, '+k')
scatter(real(zoa*p.P(p.k)/p.N), imag(zoa*p.P(p.k)/p.N), 150, 'xb');

eqpts = OA_fixedpointiteration(zoa', p);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 150, '+r');

eqpts = NewtonRaphsonIterationEntangled(zoa, p);
eqpts = eqpts(1:p.Mk) + 1i*eqpts(p.Mk+1:end);
ZOA = eqpts'*p.P(p.k)/p.N;
% scatter(real(ZOA), imag(ZOA), 150, 'xb');

% Normal OA integration
[TOA, b] = ode45(@(t,z) MFROA(t,z,p), [0, 10], zoa, opts);
Z = b*p.P(p.k)/p.N;
scatter(real(Z(1)), imag(Z(1)), 100, 'ob');
Zplot = plot(real(Z), imag(Z), 'LineWidth', 2);

% 2D OA integration
zoa2D = zeros(2*p.Mk, 1);
zoa2D(1:2:end) = real(zoa);
zoa2D(2:2:end) = imag(zoa);

[TOA, b] = ode45(@(t,z) MFROA2D(t,z,p), [0, 10], zoa2D, opts);
Z2D = b(:, 1:2:end)*p.P(p.k)/p.N + 1i*b(:, 2:2:end)*p.P(p.k)/p.N;
scatter(real(Z2D(1)), imag(Z2D(1)), 100, 'ob');
Z2Dplot = plot(real(Z2D), imag(Z2D), 'LineWidth', 2);

legend([Zplot, Z2Dplot], "OA", "OA2D", "Location", "northeast")

%%
function dfdz = MFROAJEntangled(z,p)
    x = z(1:2:end); y = z(2:2:end);
    dfdz = zeros(2*p.Mk, 2*p.Mk);
    
    etaH2k = p.eta0 + p.OA*(1 + (x.*x)/3 - 4.*x/3);
        
    for r = 1:p.Mk
        xr = x(r); yr = y(r);
        for c = 1:p.Mk
            dH2k = p.OA(r,c)*(x(c)-2)*2/3;
            if r == c
                dfdz(r + 0,c + 0) =   yr - (xr + 1)*p.delta - yr*etaH2k(r) - (xr + 1)*yr*dH2k;
                dfdz(r + 0,c + 1) =  (xr - 1) + yr*p.delta - (xr + 1)*etaH2k(r);
                dfdz(r + 1,c + 0) = -(xr - 1) - yr*p.delta + (xr + 1)*etaH2k(r) + 0.5*((xr+1)^2 - yr^2)*dH2k;
                dfdz(r + 1,c + 1) =   yr - (xr + 1)*p.delta - yr*etaH2k(r);
            else 
                dfdz(r + 0,c + 0) = -(xr + 1)*yr*dH2k;
                % dfdz(r + p.Mk,c + 0)    = 0;
                dfdz(r + 1,c + 0) = 0.5*((xr+1)^2 - yr^2)*dH2k;
                % dfdz(r + p.Mk,c + p.Mk) = 0;
            end
        end
    end
end


function [z, zs] = NewtonRaphsonIterationEntangled(z0, p)

    f = @(z, p) MFROA2D(0,z,p);    
    df = @(z, p) MFROAJEntangled(z,p);
    
    z = zeros(2*p.Mk, 1);
    z(1:2:end) = real(z0);
    z(2:2:end) = imag(z0);
    
    zc = z0';
    zval = zeros(2*p.Mk, 1);
        
    maxevals = 30;
    zs = NaN(2, maxevals+1);
    zs(1, 1) = z(1:2:end)'*p.P(p.k)/p.N;
    zs(2, 1) = z(2:2:end)'*p.P(p.k)/p.N;
    for evaltime = 1:maxevals
        z0 = z;
        
        fval = MFROA(0,zc',p);    
        zval(1:2:end) = real(fval);
        zval(2:2:end) = imag(fval);
        
        fdiv = df(z,p)\zval;
        z = z - fdiv;
        
        error = norm(z - z0);
        if error < 1.0e-9
            break
        end
        zs(1, evaltime+1) = z(1:2:end)'*p.P(p.k)/p.N;
        zs(2, evaltime+1) = z(2:2:end)'*p.P(p.k)/p.N;
        if abs(zs(1, evaltime+1) + 1i*zs(2, evaltime+1)) > 1
            warning("Out of the complex circle");
            break
        end
        
        zc = z(1:2:end) + 1i*z(2:2:end);
    end
    plot(zs(1,:), zs(2,:), 'LineWidth', 2)
    disp(['Algorithm took ', num2str(evaltime), ' steps'])
%     test = df(z0, p);
%     det(test)
%     trace(test)
%     eig(test)
end


function dfdz = MFROAJ(z,p)
    x = z(1:p.Mk); y = z(p.Mk+1:end);
    dfdz = zeros(2*p.Mk, 2*p.Mk);
    
    etaH2k = p.eta0 + p.OA*(1 + (x.*x)/3 - 4.*x/3);
        
    for r = 1:p.Mk
        xr = x(r); yr = y(r);
        for c = 1:p.Mk
            dH2k = p.OA(r,c)*(x(c)-2)*2/3;
            if r == c
                dfdz(r + 0,      c + 0) =   yr - (xr + 1)*p.delta - yr*etaH2k(r) - (xr + 1)*yr*dH2k;
                dfdz(r + 0,   c + p.Mk) =  (xr - 1) + yr*p.delta - (xr + 1)*etaH2k(r);
                dfdz(r + p.Mk,   c + 0) = -(xr - 1) - yr*p.delta + (xr + 1)*etaH2k(r) + 0.5*((xr+1)^2 - yr^2)*dH2k;
                dfdz(r + p.Mk,c + p.Mk) =   yr - (xr + 1)*p.delta - yr*etaH2k(r);
            else 
                dfdz(r + 0,   c + 0)      = -(xr + 1)*yr*dH2k;
                % dfdz(r + p.Mk,c + 0)    = 0;
                dfdz(r + p.Mk,c + 0)      = 0.5*((xr+1)^2 - yr^2)*dH2k;
                % dfdz(r + p.Mk,c + p.Mk) = 0;
            
            end
            
        end
    end
end


function [z, zs] = NewtonRaphsonIteration(z0, p)

    f = @(z, p) MFROA2D(0,z,p);    
    df = @(z, p) MFROAJ(z,p);
    
    z = zeros(2*p.Mk, 1);
    z(1:p.Mk) = real(z0);
    z(p.Mk+1:end) = imag(z0);

    maxevals = 30;
    
    zs = NaN(2, maxevals+1);
    zs(1, 1) = z(1:p.Mk)'*p.P(p.k)/p.N;
    zs(2, 1) = z(p.Mk+1:end)'*p.P(p.k)/p.N;
    for evaltime = 1:maxevals
        z0 = z;
        
%         fval = f(z,p);
%         fval = [fval(1,:), fval(2,:)]';
%         fval = df(z,p)\fval;
%         fval = [fval(1:p.Mk), fval(p.Mk+1:end)];
    
%         det(df(z,p))
        fdiv = df(z,p)\f(z,p);
        z = z - fdiv;
        
        error = norm(z - z0);
        if error < 1.0e-9
            break
        end
        zs(1, evaltime+1) = z(1:p.Mk)'*p.P(p.k)/p.N;
        zs(2, evaltime+1) = z(p.Mk+1:end)'*p.P(p.k)/p.N;
        if abs(zs(1, evaltime+1) + 1i*zs(2, evaltime+1)) > 1
            warning("Out of the complex circle");
            break
        end
    end
    plot(zs(1,:), zs(2,:), 'LineWidth', 2)
    disp(['Algorithm took ', num2str(evaltime), ' steps'])
    test = df(z0, p);
    det(test)
    trace(test)
%     eig(test)
end