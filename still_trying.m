% n = 6;
% T = diag(-4*ones(n,1)) + diag(ones(n-1,1),1)+ diag(ones(n-1,1),-1);
% A = kron(eye(n),T) + diag(ones((n)^2  - n,1),n) + diag(ones((n)^2 -n,1),-(n));

u_1 = @(x,y) x.*(x-1).*y.*(y-1);
d2u_1 = @(x,y) (2.*x.^2) + (2.*y.^2) - (2.*x) - (2.*y);
u_2 = @(x,y) sin(2*pi*x).*sin(2*pi*y);
d2u_2 = @(x,y) -8*pi^2*sin(2*pi*x).*sin(2*pi*y);
u_3 = @(x,y) (x-.5).^2 + (y-.5).^2;
d2u_3 = @(x,y) 4.*ones(size(x));
u_4 = @(x,y) (x-.5).^4 + (y-.5).^4;
d2u_4 = @(x,y) 12*(x-.5).^2 + 12*(y-.5).^2;

%% Run iterative methods
c = 1;
for i = 5%8:4:64
    n = i;
    h = 1/(n+1);
    x = 0:h:1;
    y = 0:h:1;
    [X,Y] = meshgrid(x,y); 
    sol = u_1(X,Y);
    sol1 = reshape(sol.', [],1);
    rhs = compute_gridpoints_fns(d2u_1,x,y);%reshape(twod_mult_ax(sol1, n+2,h).', [(n+2), (n+2)])';
%     [u_j,k_j] = jacobi_solve(n,rhs,sol(1,:),sol(end,:),sol(:,1),sol(:,end));
%     [u_gs,k_gs] = gauss_seidel_solve(n,rhs,sol(1,:),sol(end,:),sol(:,1),sol(:,end));
%     
    rhs2 = twod_mult_ax(sol1, n+2,h);%reshape(rhs.', [],1);%
    [u_cg, k_cg] = conjugate_gradient_test(n+2,h,reshape(rhs.', [],1),sol(1,:),sol(end,:),sol(:,1),sol(:,end));
    norm(reshape(rhs.', [],1) - rhs2,2)
% 
% % record error due to grid size
%     e_j(c) = (1/n)*norm(sol - u_j,1);
%     e_gs(c) = (1/n)*norm(sol - u_gs,1);
%     e_cg(c) = (1/n)*norm(sol - u_cg,1);
% % record iterations needed with fixed grid size
%     iter_j(c) = k_j;
%     iter_gs(c) = k_gs;
%     iter_cg(c) = k_cg;
%     c = c+1;
end

figure()
surf(X,Y,sol)
title("True Solution to -\Deltau = f");ylabel("y");xlabel("x");zlabel("u(x,y)");
ax = gca; % current axes
ax.FontSize = 14;
grid on

figure()
surf(X,Y,u_cg)
title("Conjugate Gradient to -\Deltau = f");ylabel("y");xlabel("x");zlabel("u(x,y)");
ax = gca; % current axes
ax.FontSize = 14;
grid on
