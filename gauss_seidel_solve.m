function x = gauss_seidel_solve(np, rhs)
% np - number of grid points
%  ub - upper bound
%  lb - lower bound
tol = 1e-6;
h = 1/(np-1);
x = zeros(np,np);
x_new = zeros(np,np);
r = zeros(np,np);
norm_rhs = norm( rhs(:),2);
%% Update x
tol = 1e-6;
while 1
    for i = 2:np-1
        for j = 2:np-1
            x_new(i,j) = (1/4)*(x_new(i-1,j)+x(i+1,j) + x(i,j-1) + x_new(i,j+1) - h^2 *rhs(i,j));
            r = rhs(i,j) - (1/h^2)*(x_new(i-1,j)+x(i+1,j) + x(i,j-1) + x_new(i,j+1)-4*x(i,j));
        end
    end
    if norm(r(:),2)/norm_rhs < tol
         break;
    end
    x = x_new;
end
end