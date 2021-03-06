close all;clear all;
%% Set up problems
u_1 = @(x,y) x.*(x-1).*y.*(y-1);
d2u_1 = @(x,y) 2*x.^2 + 2*y.^2 - 2*x - 2*y;
u_2 = @(x,y) sin(2*pi.*x).*sin(2*pi.*y);
d2u_2 = @(x,y) -8*pi^2*sin(2*pi.*x).*sin(2*pi.*y);
u_3 = @(x,y) (x-.5).^2 + (y-.5).^2;
d2u_3 = @(x,y) 4.*ones(size(x));
u_4 = @(x,y) (x-.5).^4 + (y-.5).^4;
d2u_4 = @(x,y) 12*(x-.5).^2 + 12*(y-.5).^2;

%% Run iterative methods
c = 1;
for i = 8:4:64
    n = i;
    h = 1/(n+1);
    x = 0:h:1;
    y = 0:h:1;
    [X,Y] = meshgrid(x,y); 
    rhs = d2u_1(X,Y);
    sol = u_1(X,Y);
%     [u_j,k_j] = jacobi_solve(n,rhs,sol(1,:),sol(end,:),sol(:,1),sol(:,end));
%     [u_gs,k_gs] = gauss_seidel_solve(n,rhs,sol(1,:),sol(end,:),sol(:,1),sol(:,end));
    
    sol1 = reshape(sol.', [],1);
    rhs1 = rhs;
    rhs1(2:n+1,2:n+1) = h^2*rhs(2:n+1,2:n+1);
    rhs1(1,:) = rhs(1,:)-sol(1,:)*h^2;
   rhs1(n+2,:) = rhs(n+2,:) - sol(n+2,:)*h^2;
   rhs1(:,1) = rhs(:,1)- sol(:,1)*h^2;
   rhs1(:,n+2) =rhs(:,n+2)- sol(:,n+2)*h^2;
    rhs1 = reshape(rhs1.', [],1);%twod_mult_ax(sol1, n+2,h);
    [u_cg, k_cg] = conjugate_gradient_test(n+2,h,rhs1,sol(1,:),sol(end,:),sol(:,1),sol(:,end));
% record error due to grid size
%     e_j(c) = (1/n)*norm(sol - u_j,1);
%     e_gs(c) = (1/n)*norm(sol - u_gs,1);
    e_cg(c) = (1/n)*norm(sol - u_cg,1);
% record iterations needed with fixed grid size
%     iter_j(c) = k_j;
%     iter_gs(c) = k_gs;
    iter_cg(c) = k_cg;
    c = c+1;
end
%% Plot convergence - number of iterations needed
% figure()
% plot(8:4:64,iter_j, 'color', 'red','LineWidth',3)
% hold on
% plot(8:4:64,iter_gs, 'o-','color', 'blue','LineWidth',3)
plot(8:4:64,iter_cg, '*-','color', 'green','LineWidth',3)
title("Number of Iterations Needed vs Grid Size- Example 2");ylabel("Iterations");xlabel("Number of Grid Points");
ax = gca; % current axes
ax.FontSize = 14;
% legend('Jacobi Method','Gauss-Seidel Method', 'Conjugate Gradient Method')
%% Plot error
% figure()
% plot(8:4:64,e_j, 'color', 'red','LineWidth',3)
% hold on
% plot(8:4:64,e_gs, 'o','color', 'blue','LineWidth',3)
%plot(8:4:64,e_cg, '*-','color', 'green','LineWidth',3)
% title("Error vs. Grid Spacing - Example 4");ylabel("Error");xlabel("Number of Grid Points");
% ax = gca; % current axes
% ax.FontSize = 14;
% legend('Jacobi Method','Gauss-Seidel Method', 'Conjugate Gradient Method')
%% Surface Plots
% figure()
% surf(X,Y,sol)
% title("True Solution to Au=\Deltau");ylabel("y");xlabel("x");zlabel("u(x,y)");
% ax = gca; % current axes
% ax.FontSize = 14;
% grid on
% 
% figure()
% surf(X,Y,u_j)
% title("Jacobi Solution to Au=\Deltau");ylabel("y");xlabel("x");zlabel("u(x,y)");
% ax = gca; % current axes
% ax.FontSize = 14;
% grid on
% 
% figure()
% surf(X,Y,u_gs)
% title("Gauss-Seidel Solution to Au=\Deltau");ylabel("y");xlabel("x");zlabel("u(x,y)");
% ax = gca; % current axes
% ax.FontSize = 14;
% grid on
% 
% figure()
% surf(X,Y,u_cg)
% title("Conjugate Gradient to Au=\Deltau");ylabel("y");xlabel("x");zlabel("u(x,y)");
% ax = gca; % current axes
% ax.FontSize = 14;
% grid on

%% Cross Section Plots
% Create x and y over the slicing plane
% xq=linspace(0,1,i);
% yq=linspace(0,1,i);
% 
% % Interpolate over the surface
% zq_j=interp2(X,Y,u_j,xq,yq); 
% zq_gs=interp2(X,Y,u_gs,xq,yq);
% zq_sol=interp2(X,Y,sol,xq,yq); 
% zq_cg=interp2(X,Y,u_cg,xq,yq); 
% figure()
% plot(xq,zq_sol,'color','black','LineWidth',3)
% hold on
% plot(xq,zq_j,'.','color','red','MarkerSize',40)
% plot(xq,zq_gs,'o','color','blue','LineWidth',3)
% plot(xq,zq_cg,'*','color', 'green','LineWidth',3)
% title("Cross Section of Solution to Au=\Deltau");xlabel("x");ylabel("u(x,y)");
% ax = gca; % current axes
% ax.FontSize = 14;
% legend('True Solution','Jacobi Method','Gauss-Seidel Method', 'Conjugate Gradient Method')