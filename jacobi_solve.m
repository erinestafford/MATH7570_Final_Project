function x = jacobi_solve(np, rhs,x0,x1,y0,y1)
% np - number of grid points
% rhs - matrix of values of 2nd derivative
% x0,x1,y0,y1 - the boundary conditions
tol = 1e-6;
h = 1/(np+1);
x = zeros(np+2,np+2);
x_new = zeros(np+2,np+2);
r = zeros(np+2,np+2);
norm_rhs = norm( rhs(:),2);
%% Update x
tol = 1e-6;
while 1
    for i = 1:np+2 %x
        for j = 1:np+2 %y
            %Apply BCs
            if(i ==1)
                x_new(i,j) = x0(j);
                r(i,j) = 0;
            elseif(i==np+2)
                x_new(i,j) = x1(j);
                r(i,j) = 0;
            elseif(j==1)
                x_new(i,j) = y0(i);
                r(i,j) = 0;
            elseif(j==np+2)
                x_new(i,j) = y1(i);
                r(i,j) = 0;
            else
                x_new(i,j) = (1/4)*(x(i-1,j)+x(i+1,j) + x(i,j-1) + x(i,j+1) - h^2 *rhs(i,j));
                r(i,j) = rhs(i,j) - (1/h^2)*(x(i-1,j)+x(i+1,j) + x(i,j-1) + x(i,j+1)-4*x(i,j));
            end
        end
    end
    if norm(r(:),2)/norm_rhs < tol
         break;
    end
    x = x_new;
end
end