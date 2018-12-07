A = [4 -1 0 4; -1 4 1 2; 0 1 4 1; 4 2 1 4];
sol = [3 1 1 2]'; %(solution)
b = A*sol
tol = 1e-3;

x = [0 0 0 0]';
r = b - A*x;
s = r;

for k = 1:100
    u = A*s;
    alpha = dot(s,r)/dot(s,u);
    x = x + alpha*s;
    r = b - A*x;
    
    if sqrt(dot(r,r))<tol
        return
    else
        beta = -dot(r,u)/dot(s,u);
        s = r + beta*s;
    end
end
