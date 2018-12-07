function x = conjugate_gradient_solve(np, rhs,x0,x1,y0,y1)
% np - number of grid points
% rhs - matrix of values of 2nd derivative
%% initialize
tol = 1e-8;
h = 1/(np+1);
x = zeros(np+2,np+2);
x_new = x; %b - Ax_0
r= h^2*rhs; %since x_0 = 0
r_new = r;
k = 0;

delta = (r(:))'*r(:);
b_delta = (r(:))'*r(:);
p = r;
%% iterate
while  norm(r(:),2)>tol
   %% calculate s_k
   s=zeros(size(p));
   for i = 2:np+1 %x
        for j = 2:np+1 %y
%             %% Apply BCs
%             if (i ==1 && j ==1) %p(1,1)
%                 s(i,j) = -4*p(i,j) + p(i,j+1) + p(i+1,j);
%             elseif (i~=1 && i~=np+2 && j==1) %column 1
%                 s(i,j) = -4*p(i,j) + p(i-1, j) + p(i,j+1) + p(i+1,j);
%             elseif (i == 1 && j ~=np+2 && j ~= 1)%row 1
%                 s(i,j) = -4*p(i,j) + p(i,j-1) + p(i,j+1) + p(i+1,j);
%             elseif (i == 1 && j==np+2 )% (1,np)
%                 s(i,j) = -4*p(i,j) + p(i+1,j) + p(i,j-1);
%             elseif ( i == np+2 && j==1 )% (np,1)
%                 s(i,j) = -4*p(i,j) + p(i-1,j) + p(i,j+1);
%             elseif ( i == np+2 && j==np+2 )% (np,np)
%                 s(i,j) = -4*p(i,j) + p(i-1,j) + p(i,j-1);
%             elseif (i ~=1 && i ~= np+2 && j==np+2) %column np
%                 s(i,j) = -4*p(i,j) + p(i-1, j) + p(i,j-1) + p(i+1,j);
%             elseif (i == np+2 && j~=np+2&&j ~=1) %row np
%                 s(i,j) = -4*p(i,j) + p(i,j-1) + p(i,j+1) + p(i-1,j);
%             else
                s(i,j) = p(i-1,j)+p(i+1,j) + p(i,j-1) + p(i,j+1)-4*p(i,j);
%             end
        end
   end
   s = s(2:np+1,2:np+1);
   alpha = dot(p(:),r(:))/dot(p(:),s(:));
   x = x + alpha*p;
   x(1,:) = x0;
   x(np+2,:) = x1;
   x(:,1) = y0;
   x(:,np+2) = y1;
   
   r = r - alpha*s;
   
   delta_new = (r(:))'*r(:);
   beta = delta_new/delta;
   p = r + beta*p;
   delta = delta_new;
   k = k+1;
end

end