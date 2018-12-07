function x = conjugate_gradient_solve(np, rhs,x0,x1,y0,y1)
% np - number of grid points
% rhs - matrix of values of 2nd derivative

%% Initialization
tol = 1e-8;
h = 1/(np+1);
x = zeros((np+2)^2,1);
r= h^2*rhs(:);
s = r;
%% Iteration

for k = 1:200
    
end


end