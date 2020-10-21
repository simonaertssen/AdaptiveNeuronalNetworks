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

%eqpts = vpasolve([f(x, y, 0, 0, 0) == 0, g(x, y, 0, 0, 0) == 0], [x,y])
disp('searching equilibria')
eqpts = solve([f(x, y, D, e, K) == 0, g(x, y, D, e, K) == 0], [x,y], 'ReturnConditions',true);
eqpts.x
eqpts.y
save('equilibria_expression', 'eqpts');
disp('equilibria found')

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



