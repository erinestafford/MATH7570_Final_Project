clear all
close all

a = 0;
b = 2*pi;
c = 1;
for N = 8:2:64
N = 8; % number of cells; N+1 is the number of grid points
L = b-a; % size of the domain
h = L/N; % grid spacing
x = linspace(a,b,N+1)'; % column vector of grid points
phi_f = @(t) t.*cos(t); % inline function for the exact solution
rho_f = @(t) 2.*sin(t) + t.*cos(t); % inline function for the exact right-hand-side

phi_exact = phi_f(x); % exact solution at the grid points

alpha = phi_f(a); % boundary condition
beta = phi_f(b); % boundary condition

phi = zeros(N-1,1); % column vector for the solution
A = zeros(N-1,N-1); % matrix of discrete Laplacian.

rho = rho_f(x(2:N)); % exact rho at the grid points
rho(1) = rho(1) + alpha/h^2; % add boundary term
rho(N-1) = rho(N-1) + beta/h^2; % add boundary term

dA = diag( 2*ones(1,N-1) ); % diagonal matrix
dAp1 = diag( -1*ones(1,N-2), 1 ); % super-diagonal matrix
dAm1 = diag( -1*ones(1,N-2), -1 ); % sub-diagonal matrix
A = (dA + dAp1 + dAm1);
A ;% display the structure of A
A = A/h^2;
%A*phi_exact'
rho;

%eigv = eig(A) % for curiosity, we check the eigenvalues of A

%% cg solve
x = zeros(length(A),1);
r = rho;
norm_rhs = norm(r,2);
d = r'*r;
bd = r'*r;
k = 0;
p = r;
tol = 1e-8;
% Iterations
while norm(r,2)/norm_rhs > tol
% s = A*p the matrix times the search direction

   s= A*p;

   a = d/(p'*s);
   x = x + a*p;
   r = r - a*s;
   if (r'*r < tol)
       break;
   end
   d_new = r'*r;
   p = r+(d_new/d)*p;
   d = d_new;
   k = k+1;
end

phi = A\rho; % solving the linear system

solver_err(c) = max( abs(A*x - rho));
c = c+1;
end
plot(8:2:64,solver_err)
%solver_err = max( abs(A*phi - rho) ) % Did the solver do his job?

%phi_num = [alpha; phi; beta]; % enlarge the solution vector by the values at the boundary

