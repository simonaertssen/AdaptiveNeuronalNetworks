%% Getting the equations
syms H(z) ODE(z, D, e, K)
H(z) = 1 + (z^2 + conj(z)^2)/6 - 4/3*real(z);
ODE(z, D, e, K) = -1i*(z - 1)^2/2 + (z + 1)^2/2*(-D + 1i*e + 1i*K*H(z));

syms f(x, y, D, e, K) g(x, y, D, e, K)
f(x, y, D, e, K) = simplify(real(ODE(x+1i*y, D, e, K)));
g(x, y, D, e, K) = simplify(imag(ODE(x+1i*y, D, e, K)));

syms J(x, y, D, e, K) [2,2]
J(x, y, D, e, K) = [diff(f(x, y, D, e, K),x), diff(f(x, y, D, e, K),y); 
                    diff(g(x, y, D, e, K),x), diff(g(x, y, D, e, K),y)];
 
syms d(x, y, D, e, K) 
d(x, y, D, e, K) = det(J(x, y, D, e, K));

syms x y D e K
assume(D,'real'); assume(e,'real'); assume(K,'real');
assume(x,'real'); assume(y,'real'); assume(x >= -1 & x <= 1); assume(y >= -1 & y <= 1);

disp('searching equilibria')
% eqpts = solve([f(x, y, D, e, K) == 0, g(x, y, D, e, K) == 0, x + y <= 1], [x,y], 'ReturnConditions',true);
% eqpts.x
% eqpts.y
% save('equilibria_expression', 'eqpts');
disp('equilibria found')

%% Testing some long functions:
syms eqptsy(y, D, e, K)
eqptsy(y, D, e, K) = 16*D^2*K^4*y^12 - 64*D*K^4*y^11 + (32*D^2*K^4 - 96*D^2*K^3*e + 96*D^2*K^3 + 64*K^4)*y^10 + (-288*D^3*K^3 - 576*D*K^4 + 288*D*K^3*e - 288*D*K^3)*y^9 + (72*D^4*K^2 + 528*D^2*K^4 + 288*D^2*K^3*e + 216*D^2*K^2*e^2 + 960*D^2*K^3 - 432*D^2*K^2*e + 216*D^2*K^2 - 192*K^3*e + 192*K^3)*y^8 + (-288*D^3*K^3 + 864*D^3*K^2*e - 1296*D^3*K^2 + 1440*D*K^3*e - 432*D*K^2*e^2 - 2208*D*K^3 + 864*D*K^2*e - 432*D*K^2)*y^7 + (1368*D^4*K^2 - 216*D^4*K*e + 512*D^2*K^4 - 1152*D^2*K^3*e - 1080*D^2*K^2*e^2 - 216*D^2*K*e^3 + 216*D^4*K + 6048*D^2*K^3 - 864*D^2*K^2*e + 648*D^2*K*e^2 + 2808*D^2*K^2 - 648*D^2*K*e + 384*K^3*e + 144*K^2*e^2 + 216*D^2*K - 288*K^2*e + 144*K^2)*y^6 + (-648*D^5*K - 4608*D^3*K^3 - 3456*D^3*K^2*e - 648*D^3*K*e^2 - 4896*D^3*K^2 + 2376*D^3*K*e - 512*D*K^4 - 4032*D*K^3*e - 864*D*K^2*e^2 + 216*D*K*e^3 - 1728*D^3*K + 1728*D*K^2*e - 648*D*K*e^2 - 2016*D*K^2 + 648*D*K*e - 216*D*K)*y^5 + (81*D^6 + 1152*D^4*K^2 + 864*D^4*K*e + 162*D^4*e^2 + 4096*D^2*K^4 + 6144*D^2*K^3*e + 3456*D^2*K^2*e^2 + 864*D^2*K*e^3 + 81*D^2*e^4 + 3672*D^4*K - 324*D^4*e + 1152*D^2*K^3 - 144*D^2*K^2*e - 1080*D^2*K*e^2 - 324*D^2*e^3 + 162*D^4 + 8064*D^2*K^2 - 1728*D^2*K*e + 486*D^2*e^2 + 1944*D^2*K - 324*D^2*e + 576*K^2*e + 81*D^2)*y^4 + (-648*D^5 - 8064*D^3*K^2 - 4752*D^3*K*e - 648*D^3*e^2 - 4752*D^3*K + 1296*D^3*e - 768*D*K^3 - 5760*D*K^2*e + 864*D*K*e^2 - 648*D^3 - 864*D*K*e)*y^3 + (864*D^4*K + 324*D^4*e + 6144*D^2*K^3 + 6912*D^2*K^2*e + 2592*D^2*K*e^2 + 324*D^2*e^3 + 1296*D^4 + 576*D^2*K^2 + 4320*D^2*K*e - 648*D^2*e^2 + 324*D^2*e)*y^2 + (432*D^3*K - 1296*D^3*e)*y - 324*D^4;

test = solve(eqptsy == 0, y, 'ReturnConditions', true)



%% Solving for the equilibria
step = 1;
Dgrid = -3:step:3;
Egrid = -10:step:10;
Kgrid = -20:step:20;

[X,Y,Z] = meshgrid(Dgrid,Egrid,Kgrid);
func = matlabFunction(f);
gunc = matlabFunction(g);

system = @(x, y, D, e, K) [func(x, y, D, e, K); gunc(x, y, D, e, K)];
% fsurf(@(x,y) func(x,y,0,0,0))

%%
% for D = Dgrid
%     for E = Egrid
%         for K = Kgrid
% %             [xval, yval] = vpasolve([f(x, y, D, E, K) == 0, g(x, y, D, E, K) == 0], [x,y]);
%             eqpt = fsolve(@(in) system(in(1), in(2), D, E, K), [0; 0], optimoptions('fsolve', 'Display', 'off'));
%             if norm(eqpt) <= 1
%             end
%         end
%     end
% end
% % eqpts = vpasolve([f(x, y, X, Y, Z) == 0, g(x, y, X, Y, Z) == 0], [x,y])

disp('done')

