function x = conjugate_gradient_solve(np, rhs,x0,x1,y0,y1)
% np - number of grid points
% rhs - matrix of values of 2nd derivative
tol = 1e-6;
h = 1/(np+1);
x = zeros(np+2,np+2);
x_new = zeros(np+2,np+2);
x_prev = zeros(np+2,np+2);
r= rhs; %since x_0 = 0
r_new = r;
r_prev = r;
while 1
    rprev_norm = norm(r_prev,2)^2;
    r_norm = norm(r,2)^2;
    for i = 1:np+2 %x
        for j = 1:np+2 %y
            %Apply BCs
            if(i ==1)
                x_new(i,j) = x0(j);
                r_new(i,j) = 0;
            elseif(i==np+2)
                x_new(i,j) = x1(j);
                r_new(i,j) = 0;
            elseif(j==1)
                x_new(i,j) = y0(i);
                r_new(i,j) = 0;
            elseif(j==np+2)
                x_new(i,j) = y1(i);
                r_new(i,j) = 0;
            else
                Ar = r_new(i-1,j)+r_new(i+1,j) + r_new(i,j-1) + r_new(i,j+1)-4*r_new(i,j);
                beta = (r(i,j)*Ar)/r_norm;
                gamma = (r_prev(i,j)*Ar)/rprev_norm;
                x_new(i,j) = x(i,j) - (r(i,j) - gamma*(x(i,j)-x_prev(i,j)))/(beta+gamma);
                r_new(i,j) = rhs(i,j) - (1/h^2)*(x(i-1,j)+x(i+1,j) + x(i,j-1) + x(i,j+1)-4*x(i,j));
            end
        end
    end
    r_prev = r;
    r = r_new;
    x_prev = x;
    x = x_new;
    if r'*r < tol
        break;
    end
end

end