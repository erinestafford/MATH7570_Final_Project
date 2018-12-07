close all;
%% Set up problems
u_1 = @(x,y) x.*(x-1).*y.*(y-1);
d2u_1 = @(x,y) 2*x.^2 + 2*y.^2 - 2*x - 2*y;
u_2 = @(x,y) sin(2*pi.*x).*sin(2*pi.*y);
d2u_2 = @(x,y) -8*pi^2*sin(2*pi.*x).*sin(2*pi.*y);
u_3 = @(x,y) (x-.5).^2 + (y-.5).^2;
d2u_3 = @(x,y) 4.*ones(size(x));
u_4 = @(x,y) (x-.5).^4 + (y-.5).^4;
d2u_4 = @(x,y) 12*(x-.5).^2 + 12*(y-.5).^2;
% c = 1;
%% Run iterative methods
for i = 55
    n = i;
    h = 1/(n+1);
    x = 0:h:1;
    y = 0:h:1;
    [X,Y] = meshgrid(x,y); 
    rhs = d2u_2(X,Y);
    sol = u_2(X,Y);
%     u_j = jacobi_solve(n,rhs,sol(1,:),sol(end,:),sol(:,1),sol(:,end));
%     u_gs = gauss_seidel_solve(n,rhs,sol(1,:),sol(end,:),sol(:,1),sol(:,end));
    u_cg = conjugate_gradient_solve(n,rhs,sol(1,:),sol(end,:),sol(:,1),sol(:,end));
%% record error due to grid size
%     e_j(c) = (1/n)*norm(sol - u_j,1);
%     e_gs(c) = (1/n)*norm(sol - u_gs,1);
%     c = c+1;
end
%% Plot error
% figure()
% plot(4:4:64,e_j)
% figure()
% plot(4:4:64,e_gs)

%% Surface Plots
% figure()
% surf(X,Y,u_j)
% title("Jacobi Solution to Au=\Deltau");ylabel("y");xlabel("x");zlabel("u(x,y)");
% ax = gca; % current axes
% ax.FontSize = 14;
% grid on
figure()
surf(X,Y,sol)
title("True Solution to Au=\Deltau");ylabel("y");xlabel("x");zlabel("u(x,y)");
ax = gca; % current axes
ax.FontSize = 14;
grid on
% figure()
% surf(X,Y,u_gs)
% title("Gauss-Seidel Solution to Au=\Deltau");ylabel("y");xlabel("x");zlabel("u(x,y)");
% ax = gca; % current axes
% ax.FontSize = 14;
% grid on

figure()
surf(X,Y,u_cg)
title("Conjugate Gradient to Au=\Deltau");ylabel("y");xlabel("x");zlabel("u(x,y)");
ax = gca; % current axes
ax.FontSize = 14;
grid on

%% Cross Section Plots
% Create x and y over the slicing plane
% xq=linspace(0,1,i);
% yq=linspace(0,1,i);

% Interpolate over the surface
% zq_j=interp2(X,Y,u_j,xq,yq); 
% zq_gs=interp2(X,Y,u_gs,xq,yq);
% zq_sol=interp2(X,Y,sol,xq,yq); 
% figure()
% plot(xq,zq_j,'LineWidth',3)
% title("Cross Section of Jacobi Solution to Au=\Deltau");xlabel("x");ylabel("u(x,y)");
% ax = gca; % current axes
% ax.FontSize = 14;
% figure()
% plot(xq,zq_gs,'LineWidth',3)
% title("Cross Section of Gauss-Seidel Solution to Au=\Deltau");xlabel("x");ylabel("u(x,y)");
% ax = gca; % current axes
% ax.FontSize = 14;
% figure()
% plot(xq,zq_sol,'LineWidth',3)
% title("Cross Section of True Solution to Au=\Deltau");xlabel("x");ylabel("u(x,y)");
% ax = gca; % current axes
% ax.FontSize = 14;
