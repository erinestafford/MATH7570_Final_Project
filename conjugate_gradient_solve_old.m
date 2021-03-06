function x = conjugate_gradient_solve_old(np, rhs,x0,x1,y0,y1)
% np - number of grid points
% rhs - matrix of values of 2nd derivative
%% initialize
tol = 1e-8;
h = 1/(np+1);
x = zeros(np+2,np+2);
% x_new = x; %b - Ax_0
r= h^2*rhs; %since x_0 = 0
% r_new = r;

k = 0;
r_d = r(2:np+1);
delta = (r_d(:))'*r_d(:);
b_delta = (r_d(:))'*r_d(:);

p = r;
%% iterate
while  k <1000
   %% calculate s_k
   s=mult_ax(p(:),np);
   s = reshape(s,[np+2,np+2]);
%    for i = 2:np+1 %x
%         for j = 2:np+1 %y
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
%                 s(i,j) = p(i-1,j)+p(i+1,j) + p(i,j-1) + p(i,j+1)-4*p(i,j);
%             end
%         end
%    end
   

   
   alpha = dot(p(:),r(:))/dot(p(:),s(:));
   
   x = x + alpha*p;
   x(1,:) = x0;
   x(np+2,:) = x1;
   x(:,1) = y0;
   x(:,np+2) = y1;
   
   ax = mult_ax(x(:),np);
   ax = reshape(ax,[np+2,np+2]);
   r = h^2*rhs - ax;
   
   if sqrt(dot(r(:),r(:)))<tol
        return
    else
        beta = -dot(r(:),s(:))/dot(p(:),s(:));
        p = r + beta*p;
    end
   k = k+1;
end

end