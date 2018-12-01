close all;
u_1 = @(x,y) x.*(x-1).*y.*(y-1);
d2u_1 = @(x,y) 2*x.^2 + 2*y.^2 - 2*x - 2*y;
u_2 = @(x,y) sin(2*pi.*x).*sin(2*pi.*y);
d2u_2 = @(x,y) -8*pi^2*sin(2*pi.*x).*sin(2*pi.*y);
c = 1;
for i = 2:2:64
    n = i;
    h = 1/(n+1);
    x = 0:h:1;
    y = 0:h:1;
    [X,Y] = meshgrid(x,y); %X and Y are (n+2)x (n+2)
    rhs = d2u_1(X,Y);
    
    u1_j = jacobi_solve(n,rhs);
    sol = u_1(X,Y);
%     u1_gs = gauss_seidel_solve(n,rhs);

    e_j(c) = sqrt((1/n)*sum(sum((sol - u1_j).^2)));
%     e_gs(c) = sqrt((1/n)*norm(u_1(X,Y) - u1_gs,2));
    c = c+1;
end
figure()
plot(2:2:64,e_j)
% figure()
% plot(4:2:64,e_gs)

% sol = u_1(X,Y);
% figure()
% plot3(X,Y,u1_j)
% figure()
% plot3(X,Y,sol)
% grid on
% figure()
% plot3(X,Y,u1_gs)
% grid on