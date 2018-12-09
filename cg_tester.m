clear all; close all;
u_1 = @(x,y) x.*(x-1).*y.*(y-1);
%d2u_1 = @(x,y) -2.*ones(size(x));

n = 4;
h = 1/(n-1);
x = 0:h:1;
y = 0:h:1;
[X,Y] = meshgrid(x,y); 
sol1 = u_1(X,Y);
sol = reshape(sol1.', [],1);
rhs = twod_mult_ax(sol, n,h);
[u_cg, k_cg] = conjugate_gradient_test(n,h,rhs,sol(1,:),sol(end,:),sol(:,1),sol(:,end));

figure()
subplot(2,2,1)
surf(X,Y,u_cg)
title("Conjugate Gradient to Au=\Deltau");ylabel("y");xlabel("x");zlabel("u(x,y)");
ax = gca; % current axes
ax.FontSize = 14;
grid on
subplot(2,2,2)
surf(X,Y,sol1)
title("True Solution to Au=\Deltau");ylabel("y");xlabel("x");zlabel("u(x,y)");
ax = gca; % current axes
ax.FontSize = 14;
grid on