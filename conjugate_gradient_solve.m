function [x,k] = conjugate_gradient_solve(np, rhs,x0,x1,y0,y1)
% np - number of grid points
% rhs - matrix of values of 2nd derivative

%% Initialization
tol = 1e-8;
h = 1/(np+1);
x = zeros(np+2,np+2); %initial guess is zeros
x_p = x;
x_n = x;
r= rhs;
r_p = r;
%% Iteration
k = 0;
while (k < (np+2)^2 && r(:)'*r(:)> tol)
    %% Ar = A*r the matrix times the search direction
  for i = 1:np+2 %x
        for j = 1:np+2 %y
            Ar(i,j) = 0;
            % Apply BCs
            if (i ==1 && j ==1) %p(1,1)
                Ar(i,j) = -4*r(i,j) + r(i,j+1) + r(i+1,j);
            elseif (i~=1 && i~=np+2 && j==1) %column 1
                Ar(i,j) = -4*r(i,j) + r(i-1, j) + r(i,j+1) + r(i+1,j);
            elseif (i == 1 && j ~=np+2 && j ~= 1)%row 1
                Ar(i,j) = -4*r(i,j) + r(i,j-1) + r(i,j+1) + r(i+1,j);
            elseif (i == 1 && j==np+2 )% (1,np)
                Ar(i,j) = -4*r(i,j) + r(i+1,j) + r(i,j-1);
            elseif ( i == np+2 && j==1 )% (np,1)
                Ar(i,j) = -4*r(i,j) + r(i-1,j) + r(i,j+1);
            elseif ( i == np+2 && j==np+2 )% (np,np)
                Ar(i,j) = -4*r(i,j) + r(i-1,j) + r(i,j-1);
            elseif (i ~=1 && i ~= np+2 && j==np+2) %column np
                Ar(i,j) = -4*r(i,j) + r(i-1, j) + r(i,j-1) + r(i+1,j);
            elseif (i == np+2 && j~=np+2&&j ~=1) %row np
                Ar(i,j) = -4*r(i,j) + r(i,j-1) +r(i,j+1) + r(i-1,j);
            else
                Ar(i,j) = (r(i-1,j)+r(i+1,j) + r(i,j-1) + r(i,j+1)-4*r(i,j));
            end
        end
   end
   Ar = Ar/h^2;
    %% Calculate beta and gamma
    beta = (r(:)'*Ar(:))/(norm(r(:),2)^2);
    gamma = (r_p(:)'*Ar(:))/(norm(r_p(:),2)^2);
    %% Update x
    x_n = x + (1/(beta+gamma))*(r - gamma*(x-x_p));
    x_n(1,:) = x0;
    x_n(np+2,:) = x1;
    x_n(:,1) = y0;
    x_n(:,np+2) = y1;
    
    x_p = x;
    x = x_n;
    
    %% get alpha
    alpha = -(beta + gamma);
    %% get new r
    new_r = (1/alpha)*(Ar - beta*r - gamma*r_p);
    new_r(1,:) = zeros(size(new_r(1,:)));
    new_r(np+2,:) = zeros(size(new_r(np+2,:)));
    new_r(:,1) = zeros(size(new_r(:,1)));
    new_r(:,np+2) = zeros(size(new_r(:,np+2)));
    
    r_p = r;
    r = new_r;
    
   %% update number of iterations
   k = k+1;
   
   
end


end