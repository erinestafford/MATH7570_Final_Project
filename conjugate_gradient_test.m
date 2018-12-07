function x = conjugate_gradient_test
% CONJGRAD  Conjugate Gradient Method.
%   X = CONJGRAD(A,B) attemps to solve the system of linear equations A*X=B
%   for X. The N-by-N coefficient matrix A must be symmetric and the right
%   hand side column vector B must have length N.
%
%   X = CONJGRAD(A,B,TOL) specifies the tolerance of the method. The
%   default is 1e-10.
%
% Example (highlight lines between %{ and %}, then press F9):
u_1 = @(x,y) x.*(x-1).*y.*(y-1);
d2u_1 = @(x,y) 2*x.^2 + 2*y.^2 - 2*x - 2*y;
n = 50;
h = 1/(n+1);
x = 0:h:1;
y = 0:h:1;
[X,Y] = meshgrid(x,y); 
b = d2u_1(X,Y);
b = h^2*b(:);
sol = u_1(X,Y);
sol = sol(:);
A = zeros((n+2)^2,(n+2)^2);
A(1:n,1:n) = diag(-4*ones(n,1)) + diag(ones(n-1,1),-1) + diag(ones(n-1,1),1);
A(n+1:2*n,n+1:2*n) = diag(-4*ones(n,1)) + diag(ones(n-1,1),-1) + diag(ones(n-1,1),1);
A(2*n+1:3*n,2*n+1:3*n) = diag(-4*ones(n,1)) + diag(ones(n-1,1),-1) + diag(ones(n-1,1),1);
A(3*n+1:4*n,3*n+1:4*n) = diag(-4*ones(n,1)) + diag(ones(n-1,1),-1) + diag(ones(n-1,1),1);
A((n+1):n^2,1:n) = eye(size(A((n+1):n^2,1:n)));
A(1:n,(n+1):n^2) = eye(size(A(1:n,(n+1):n^2)));

b-A*sol
    if nargin<3
        tol=1e-10;
    end
    x = b;
    r = b - A*x;
    if norm(r) < tol
        return
    end
    y = -r;
    z = A*y;
    s = y'*z;
    t = (r'*y)/s;
    x = x + t*y;
  
    for k = 1:numel(b);
       r = r - t*z;
       if( norm(r) < tol )
            return;
       end
       B = (r'*z)/s;
       y = -r + B*y;
       z = A*y;
       s = y'*z;
       t = (r'*y)/s;
       x = x + t*y;
    end
 end
