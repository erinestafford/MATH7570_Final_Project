n = 5;
I = eye(n-1);
T = diag(-4*ones(n-1,1)) + diag(ones(n-2,1),-1)+ diag(ones(n-2,1),1);
A = [T  I zeros(size(I)) zeros(size(I)); I T I zeros(size(I)); zeros(size(I)) I T I;zeros(size(I)) zeros(size(I)) I T];

u_1 = @(x,y) x*(x-1)*y*(y-1);
u_2 = @(x,y) sin(2*pi*x)*sin(2*pi*y);
