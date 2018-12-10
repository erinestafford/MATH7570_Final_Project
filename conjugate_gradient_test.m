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
bd = rhs'*rhs;
k = 0;
p = r;

%% Iterations
while norm(r,2)/norm_rhs > tol
%% s = A*p the matrix times the search direction

   s= twod_mult_ax(p, np,h);

   a = d/(p'*s);
   x = x + a*p;
   %% Apply BCs 
   x = reshape(x.', [np, np])';
   x(1,:) = rhs_mat(1,:)-x0/h^2;
   x(np,:) = rhs_mat(np,:) - x1/h^2;
   x(:,1) = rhs_mat(:,1)- y0/h^2;
   x(:,np) =rhs_mat(:,np)- y1/h^2;
   
   x = reshape(x.', [],1);
   %%
   r = r - a*s;
   if (r'*r < tol)
       break;
   end
   d_new = r'*r;
   p = r+(d_new/d)*p;
   d = d_new;
   k = k+1;
end
x = reshape(x.', [np,np]);
end