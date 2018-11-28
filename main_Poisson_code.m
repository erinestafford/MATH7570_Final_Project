n = 9;
u_1 = @(x,y) x.*(x-1).*y.*(y-1);
d2u_1 = @(x,y) 2*x.^2 + 2*y.^2 - 2*x - 2*y;
u_2 = @(x,y) sin(2*pi.*x).*sin(2*pi.*y);

x = 0:1/(n-1):1;
y = 0:1/(n-1):1;
[X,Y] = meshgrid(x,y);
rhs = d2u_1(X,Y);
jacobi_solve(n,rhs);
gauss_seidel_solve(n,rhs)
u_1(X,Y)