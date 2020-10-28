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
pars.N = 1000;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
pars.eta0 = 0.5; pars.delta = 0.7; pars.K = 2;
pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;

seed = 2; rng(seed);
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
netp = 0.1;
p = prepareOAparameters(make_randomparameters(pars, 0.44));

%% Test how we can get the fixpoints of the system: this works!
clc

figure; hold on; box on; grid on; axis square;
phasespaceplot();
drawdiraclimitcycle();

% One way:
% IC = -0.8*ones(p.N,1);
% OAIC = zeros(1,p.l);
% for i = 1:p.l
%     OAIC(i) = sum(exp(1i*IC(p.degrees_i == p.k(i)))) / (p.P(p.k(i))+1.0e-24);
% end
% scatter(real(orderparameter(IC)), imag(orderparameter(IC)), 150, '+')
% scatter(real(OAIC*p.P(p.k)/p.N), imag(OAIC*p.P(p.k)/p.N), 150, 'x')
% 
% eqpts = findeqpts(OAIC', p);
% ZOA = eqpts'*p.P(p.k)/p.N;
% scatter(real(ZOA), imag(ZOA), 150, 'xr');

% The other way:
z0 = -0.9 - 1i*0.2;
OAIC = find_ICs(z0*ones(1, p.l), z0, p.P(p.k)/p.N);
scatter(real(z0), imag(z0), 150, '+k')
scatter(real(OAIC*p.P(p.k)/p.N), imag(OAIC*p.P(p.k)/p.N), 150, 'xg');
eqpts = findeqpts(OAIC', p);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 150, 'xr');
[~, Z] = OA_simulatenetwork(0, 7, z0, p);
scatter(real(Z(1)), imag(Z(1)), 100, 'ob', 'filled')
plot(real(Z), imag(Z))

z0 = -0.6 - 1i*0.4;
IC = find_ICs(z0*ones(1, p.l), z0, p.P(p.k)/p.N);
scatter(real(z0), imag(z0), 150, '+')
scatter(real(IC*p.P(p.k)/p.N), imag(IC*p.P(p.k)/p.N), 150, 'x');
eqpts = findeqpts(IC', p);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 150, 'xr');
[~, Z] = OA_simulatenetwork(0, 7, z0, p);
scatter(real(Z(1)), imag(Z(1)), 100, 'ob', 'filled')
plot(real(Z), imag(Z))

z0 = -0.1 - 1i*0.1;
IC = find_ICs(z0*ones(1, p.l), z0, p.P(p.k)/p.N);
scatter(real(z0), imag(z0), 150, '+')
scatter(real(IC*p.P(p.k)/p.N), imag(IC*p.P(p.k)/p.N), 150, 'x');
eqpts = findeqpts(IC', p);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 150, 'xr');
[~, Z] = OA_simulatenetwork(0, 7, z0, p);
scatter(real(Z(1)), imag(Z(1)), 100, 'ob', 'filled')
plot(real(Z), imag(Z))

%% Try with fsolve?

figure; hold on; box on; grid on; axis square;
phasespaceplot();
drawdiraclimitcycle();

z0 = 0.8 - 1i*0.5;
OAIC = find_ICs(z0*ones(1, p.l), z0, p.P(p.k)/p.N);
scatter(real(z0), imag(z0), 150, '+k')
scatter(real(OAIC*p.P(p.k)/p.N), imag(OAIC*p.P(p.k)/p.N), 150, 'x');

eqpts = findeqpts(OAIC', p);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 150, 'xr');

opts = optimset('Display','off', 'Algorithm', 'levenberg-marquardt');
[eqpts,~,~,~,JACOB] = fsolve(@(x) MFROA(0,x,p), OAIC', opts);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 150, 'or');

[~, Z] = OA_simulatenetwork(0, 7, z0, p);
scatter(real(Z(1)), imag(Z(1)), 100, 'ob', 'filled')
plot(real(Z), imag(Z))

%% Make a Newton-Raphson iteration
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

clc; hold on

z0 = 0.8 - 1i*0.5;
OAIC = find_ICs(z0*ones(1, p.l), z0, p.P(p.k)/p.N);
scatter(real(z0), imag(z0), 150, '+k')
scatter(real(OAIC*p.P(p.k)/p.N), imag(OAIC*p.P(p.k)/p.N), 150, 'x');

eqpts = NewtonRaphsonIteration(OAIC', p);
ZOA = eqpts'*p.P(p.k)/p.N
scatter(real(ZOA), imag(ZOA), 150, 'xr');

[~, Z] = OA_simulatenetwork(0, 7, z0, p);
scatter(real(Z(1)), imag(Z(1)), 100, 'ob')
plot(real(Z), imag(Z))

phasespaceplot();


function x = NewtonRaphsonIteration(x0, p)
    xs = x0'*p.P(p.k)/p.N;
    
    f = @(z, p) MFROA(0,z,p);
    function dfdz = df(z, p)
        zc = conj(z);
        H = (1 + (z.*z + zc.*zc)/6 - 4.*real(z)/3);
        I = -p.delta + 1i*p.eta0 + 1i*p.K*p.OA*H;
        dIdz = 1i*p.K*p.OA*(z - 2)/3;
        dfdz = -1i*(z-1) + (z+1).*I + 1/2*(z+1).*(z+1).*dIdz;
    end
    size(f(x0,p))
    size(df(x0,p))
    size(f(x0,p)'/df(x0,p))
    return
    x = x0;
    maxevals = 400;
    for evaltime = 1:maxevals
        x0 = x;
        x = x - (f(x,p)'/df(x,p))';
        size(x)
%         if norm(x) > 1
%             x = 1./x;
%         end
        error = norm(x - x0);
        if error < 1.0e-7
            break
        end
        xs = [xs, x'*p.P(p.k)/p.N];
        plot(real(xs), imag(xs))
    end
    disp(['Algorithm took ', num2str(evaltime), ' steps'])
end

%% Functions
function bhat = findeqpts(b0, p)
    bs = b0'*p.P(p.k)/p.N;
    ntimes = 0;
    bhat = b0; bhatold = -b0;

    while norm(bhat - bhatold) > 1.0e-22 && ntimes < 400
        ntimes = ntimes + 1;
        bhatold = bhat;
        
        bhatc = conj(bhat);
        H = (1 + (bhat.*bhat + bhatc.*bhatc)/6 - 4.*real(bhat)/3);
        z = sqrt(-1i*p.delta +    p.eta0 +    p.K*p.OA*H);
        z = sqrt(   (p.delta + 1i*p.eta0 + 1i*p.K*p.OA*H)/1i);
        bhat = (1 + z)./(1 - z);
        
        if norm(bhat) > 1
            bhat = (1 - z)./(1 + z);
        end
        bs = [bs, bhat'*p.P(p.k)/p.N];
        
    end
    disp(['Algorithm took ', num2str(ntimes), ' steps'])
    plot(real(bs), imag(bs))
end

function bhat = findeqptsreversed(b0, p)
    bs = b0'*p.P(p.k)/p.N;
    ntimes = 0;
    bhat = b0; bhatold = -b0;

    while norm(bhat - bhatold) > 1.0e-22 && ntimes < 400
        ntimes = ntimes + 1;
        bhatold = bhat;
        
        z = (1 + bhat)./(1 - bhat);
        if norm(z) < 1
            bhat = (1 - bhat)./(1 + bhat);
        end
        
        zc = conj(z);
        H = (1 + (z.*z + zc.*zc)/6 - 4.*real(z)/3);
        bhat = sqrt(-1i*p.delta +    p.eta0 +    p.K*p.OA*H);        
        
        bs = [bs, bhat'*p.P(p.k)/p.N];
        
    end
    disp(['Algorithm took ', num2str(ntimes), ' steps'])
    plot(real(bs), imag(bs))
end