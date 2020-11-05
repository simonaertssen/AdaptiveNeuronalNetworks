close all; clear all; clc;
addpath('../Functions');
addpath('../Mean Field Reductions/');

%% Test fixpoint iteration of the logistic map
logistic = @(x) 2.8*x.*(1-x);

xold = 0; xnew = 0.1;
times = 0;
while norm(xnew - xold) > 1.0e-12 && times < 100
    times = times + 1;
    xold = xnew;
    xnew = logistic(xnew);
end

hold on
x = 0:0.01:1;
plot(x, logistic(x))
plot(x, x)
scatter(xnew, logistic(xnew))

%% Parameters
pars.N = 5000;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
pars.eta0 = 0.4; pars.delta = 0.7; pars.K = 2;
% pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;

seed = 1; rng(seed);
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
p = prepareOAparameters(make_randomparameters(pars, 0.33));

%% Test whether we can see the phase space of the system with vectors:
figure; hold on; box on; grid on; axis square;

startx = [0, -0.8, -0.6, 0, 0]; starty = [0, 0.2, 0.4, -1, -0.74];
tlengths = [1.5, 0.5, 0.65, 2.6, 2.2];
bw = -0.5;
opts = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);

for i = 1:length(startx)
    IC = startx(i) + starty(i)*1i;
    zoa = ones(p.Mk,1)*IC;
    zoa = find_ICs(IC*ones(1, p.Mk), IC, p.P(p.k)/p.N);
    [~, b_s] = ode45(@(t,x) MFROA(t,x,p), [0 bw], zoa, opts);
    [t, b_s] = ode45(@(t,x) MFROA(t,x,p), [bw, tlengths(i)], b_s(end,:), opts);
    transients = find(t >= 0, 1, 'first');
    Z = b_s * p.P(p.k)/p.N;
%     [TOA, Z] = OA_simulatenetwork(0, tlengths(i), OAIC, p);
%     transients = 1;
    plot(real(Z(transients:end)), imag(Z(transients:end)), 'LineWidth', 2);
    scatter(real(Z(transients)), imag(Z(transients)), 50, [0 0.3070 0.5010], 'filled', 'o', 'LineWidth',2);
    disp('done')
end
phasespaceplot();

%% Test how we can get the fixpoints of the system: this works!
clc

figure; hold on; box on; grid on; axis square;
phasespaceplot();
drawfixeddegreelimitcycle();

% The other way:
Z0 = 0.7 - 1i*0.2;
zoa = map_Ztozoa(Z0,p)';
scatter(real(Z0), imag(Z0), 150, '+k')
scatter(real(zoa*p.P(p.k)/p.N), imag(zoa*p.P(p.k)/p.N), 150, 'xg');
eqpts = findeqpts(zoa', p);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 150, 'xr');
[~, Z] = OA_simulatenetwork(0, 7, zoa, p);
scatter(real(Z(1)), imag(Z(1)), 100, 'ob', 'filled')
plot(real(Z), imag(Z))

Z0 = -0.6 - 1i*0.4;
zoa = map_Ztozoa(Z0,p)';
scatter(real(Z0), imag(Z0), 150, '+')
scatter(real(zoa*p.P(p.k)/p.N), imag(zoa*p.P(p.k)/p.N), 150, 'x');
eqpts = findeqpts(zoa', p);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 150, 'xr');
[~, Z] = OA_simulatenetwork(0, 7, zoa, p);
scatter(real(Z(1)), imag(Z(1)), 100, 'ob', 'filled')
plot(real(Z), imag(Z))

Z0 = -0.1 - 1i*0.1;
zoa = map_Ztozoa(Z0,p)';
scatter(real(Z0), imag(Z0), 150, '+')
scatter(real(zoa*p.P(p.k)/p.N), imag(zoa*p.P(p.k)/p.N), 150, 'x');
eqpts = findeqpts(zoa', p);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 150, 'xr');
[~, Z] = OA_simulatenetwork(0, 7, zoa, p);
scatter(real(Z(1)), imag(Z(1)), 100, 'ob', 'filled')
plot(real(Z), imag(Z))

%% Try with fsolve?

figure; hold on; box on; grid on; axis square;
phasespaceplot();
drawfixeddegreelimitcycle();

Z0 = 0.8 - 1i*0.5;
zoa = map_Ztozoa(Z0, p)';
scatter(real(Z0), imag(Z0), 150, '+k')
scatter(real(zoa*p.P(p.k)/p.N), imag(zoa*p.P(p.k)/p.N), 150, 'x');

eqpts = findeqpts(zoa', p);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 150, 'xr');

opts = optimset('Display','off', 'Algorithm', 'levenberg-marquardt');
% [eqpts,~,~,~,JACOB] = fsolve(@(x) MFROA2D(0,x,p), [real(zoa'), imag(zoa')], opts);
size(Z0)
size(map_Ztozoa(Z0, p))
size(MFROA(0,map_Ztozoa(Z0, p),p))
size(map_zoatoZ(MFROA(0,map_Ztozoa(Z0, p),p)',p))

[ZOA,~,~,~,JACOB] = fsolve(@(x) map_zoatoZ(MFROA(0,map_Ztozoa(Z0, p),p)',p), Z0, opts);

% ZOA = eqpts*p.P(p.k)/p.N
scatter(real(ZOA), imag(ZOA), 200, 'xk', 'LineWidth',2);

[~, Z] = OA_simulatenetwork(0, 7, zoa, p);
scatter(real(Z(1)), imag(Z(1)), 100, 'ob', 'filled')
plot(real(Z), imag(Z))

%% Make a Newton-Raphson iteration on our system
% The function is the following:
% dz/dt = -1i/2*(z-1)^2 + 1/2*(z+1)^2*I
% I = -delta + 1i*eta0 + 1i*K*OA*H/meandegree)
% H = [1 + 1/6*(z^2+zc^2) - 4/3*real(z)]

% Remember that if z = x + 1i*y then zc = x - 1i*y.
% Re(z) = (z + zc)/2 and Im(z) = -1i*(z - zc)/2.
% partd/partdz = partdx/partdz*partd/partdx + partdy/partdz*partd/partdy
% with partdx/partdz = 1/2 and partdy/partdz = -1i/2:
% partd/partdz = 1/2*(partd/partdx - 1i*partd/partdy) -> Wirtinger operators
% So partdzc/partdz = 1/2*(partd(zc)/partdx - 1i*partd(zc)/partdy)
%                   = 1/2*(1 - 1i*(-1i)) = 0
% See https://www.researchgate.net/publication/317953764_Derivative_of_Complex_Conjugate_and_Magnitude_-_Rev_1

% zc^2 = (x - 1i*y)*(x - 1i*y) = x^2 - 2*1i*x*y + (-1)*y^2 = x^2 - y^2 - 2*1i*x*y
% Then partd(zc^2)/partdz = 1/2*(partd(zc^2)/partdx - 1i*partd(zc^2)/partdy)
%                         = 1/2*(2*x - 2*1i*y - 1i*(-2*y - 2*1i*x))
%                         = x - 1i*y - 1i*(-y - 1i*x) = x - 1i*y + 1i*y + (-1)*x
%                         = 0

% Then partd(Re(z))/partdz = 1/2*(partd(x)/partdx - 1i*partd(x)/partdy)
%                          = 1/2*(1)
%                          = 1/2


% The derivative is now found as:
% d(dz/dt)/dz = -1i/2*2*(z-1) + 1/2*2*(z+1)*I + 1/2*(z+1)^2*dI/dz
%             = -1i*(z-1) + (z+1)*I + 1/2*(z+1)^2*dI/dz
% dI/dz = 1i*K*OA*(dH/dz)/meandegree
% dH/dz = 2/6*(z) - 4/3*1/2 = 1/3*z - 2/3

figure; hold on
phasespaceplot();

Z0 = 0.3 - 1i*0.5;
zoa = map_Ztozoa(Z0, p)';
scatter(real(Z0), imag(Z0), 150, '+k')
scatter(real(zoa*p.P(p.k)/p.N), imag(zoa*p.P(p.k)/p.N), 150, 'xb');

eqpts = findeqpts(zoa', p);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 150, '+r');

eqpts = NewtonRaphsonIteration(zoa', p);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 150, 'xb');

[~, Z] = OA_simulatenetwork(0, 10, zoa, p);
scatter(real(Z(1)), imag(Z(1)), 100, 'ob');
plot(real(Z), imag(Z));



%% Functions:

function dzdt = MFROA2D(t, z, p)
    Z = MFROA(t, z, p);
    dzdt = [real(Z), imag(Z)];
end

function [x, xs] = NewtonRaphsonIteration(x0, p)

    f = @(z, p) MFROA(0,z,p);
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


%% Try a Newton-Raphson system iteration:
% testNR([0.83;0.03])

function x = testNR(x0)
    function out = f(x)
%         out = zeros(3,1);
%         out(1) = 3*x(1) - cos(x(2)*x(3)) - 3/2;
%         out(2) = 4*x(1)*x(1) - 625*x(2)*x(2) + 2*x(3) - 1;
%         out(3) = x(3) + exp(-x(1)*x(2)) + 9;
        out = zeros(2,1);
        out(1) = -x(1)^3 + x(2);
        out(2) = x(1)^2 + x(2)^2-1;
    end

    function out = J(x)
%         out = zeros(3,3);
%         out(1,1) = 3; out(1,2) = x(3)*sin(x(2)*x(3)); out(1,3) = x(2)*sin(x(2)*x(3));
%         out(2,1) = 8*x(1); out(2,2) = -1250*x(2); out(1,3) = 2;
%         out(3,1) = -x(2)*exp(-x(1)*x(2)); out(3,2) = -x(1)*exp(-x(1)*x(2)); out(3,3) = 20;
%         out = [3, x(3)*sin(x(2)*x(3)), x(2)*sin(x(2)*x(3));
%                8*x(1), -1250*x(2), 2;
%                -x(2)*exp(-x(1)*x(2)), -x(1)*exp(-x(1)*x(2)), 20];
        out(1,1) = -3*x(1)^2; 
        out(1,2) = 1;
        out(2,1) = 2*x(1);
        out(2,2) = 2*x(2);
    end
    
    x = x0;
    maxevals = 1000;
    for evaltime = 1:maxevals
        x0 = x;
       
        x = x - J(x)\f(x);
        error = norm(x - x0);
        if error < 1.0e-24
            break
        end
    end
end

function bhat = findeqpts(b0, p)
    bs = b0'*p.P(p.k)/p.N;
    ntimes = 0;
    bhat = b0; bhatold = -b0;

    while norm(bhat - bhatold) > 1.0e-22 && ntimes < 400
        ntimes = ntimes + 1;
        bhatold = bhat;
        
        bhatc = conj(bhat);
        H = (1 + (bhat.*bhat + bhatc.*bhatc)/6 - 4.*real(bhat)/3);
        z = sqrt(-1i*p.delta +    p.eta0 +    p.OA*H);
        z = sqrt(   (p.delta + 1i*p.eta0 + 1i*p.OA*H)/1i);
        bhat = (1 + z)./(1 - z);
        
        if norm(bhat) > 1
            bhat = (1 - z)./(1 + z);
        end
        bs = [bs, bhat'*p.P(p.k)/p.N];
        
    end
    disp(['Algorithm took ', num2str(ntimes), ' steps'])
    plot(real(bs), imag(bs))
end

