function [x,k] = conjugate_gradient_test(np,h, rhs,x0,x1,y0,y1)
% np - number of grid points
% rhs - matrix of values of 2nd derivative

%% Initialization
tol = 1e-8;
x = zeros((np)^2,1); %initial guess is zeros
r= rhs;
rhs_mat = reshape(x.', [np, np])';
norm_rhs = norm(r,2);
d = r'*r;
dp = d;
bd = rhs'*rhs;
k = 0;
p = r;
n = np;
T = diag(-4*ones(n,1)) + diag(ones(n-1,1),1)+ diag(ones(n-1,1),-1);
A = kron(eye(n),T) + diag(ones((n)^2  - n,1),n) + diag(ones((n)^2 -n,1),-(n));

while sqrt(d) > tol*norm(rhs,2) && k <=(np^2)
    k = k+1;
    w = 1/h^2*A*p;
    a = dot(r,r)/dot(p,w);
    x = x + a*p;
    x = reshape(x, [np, np])';
   x(1,:) = x0;
   x(np,:) =x1;
   x(:,1) = y0;
   x(:,np) =y1;
   
   x = reshape(x, [],1);
    r = (rhs-A*x) - a*w;
    p  = r+(-dot(r,w)/dot(p,w))*p;
end
x = reshape(x, [np,np]);
end