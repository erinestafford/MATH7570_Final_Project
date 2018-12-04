close all;
u_1 = @(x,y) x.*(x-1).*y.*(y-1);
d2u_1 = @(x,y) 2*x.^2 + 2*y.^2 - 2*x - 2*y;
u_2 = @(x,y) sin(2*pi.*x).*sin(2*pi.*y);
d2u_2 = @(x,y) -8*pi^2*sin(2*pi.*x).*sin(2*pi.*y);
u_3 = @(x,y) (x-.5).^2 + (y-.5).^2;
d2u_3 = @(x,y) 4.*ones(size(x));
c = 1;
for i = 55
    n = i;
    h = 1/(n+1);
    x = 0:h:1;
    y = 0:h:1;
    [X,Y] = meshgrid(x,y); 
    rhs = d2u_3(X,Y);
    
    sol = u_3(X,Y);
    u_j = jacobi_solve(n,rhs,sol(1,:),sol(end,:),sol(:,1),sol(:,end));
    u_gs = gauss_seidel_solve(n,rhs,sol(1,:),sol(end,:),sol(:,1),sol(:,end));

    e_j(c) = (1/n)*norm(sol - u_j,1);
%     e_gs(c) = (1/n)*norm(sol - u_gs,1);
    c = c+1;
end
figure()
plot(4:4:64,e_j)
% figure()
% plot(4:4:64,e_gs)

figure()
plot3(X,Y,u_j)
grid on
figure()
plot3(X,Y,sol)
grid on
% figure()
% plot3(X,Y,u1_gs)
% grid on