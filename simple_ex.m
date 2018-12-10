%spy(A)
A = [5 2 1; 2 10 2; 1 2 15];
x = ones(3,1);
b = A * x;
r = b;
tol = 1e-10;
k = 1;
x =zeros(3,1); 
delta = r'*r;
b_d = b'*b;
p = r;
% while r'*r > tol
for i = 1:20
    s = A*p
    p
    a = delta/(p'*s)
    x = x+a*p
    r = r-a*s
    if( r'*r < tol)
        break;
    end
    delta_n = r'*r;
    p = r + (delta_n/delta)*p;
    delta = delta_n;
    
end
x
% n = 4;
% T = diag(-4*ones(4,1)) + diag(ones(3,1),1)+ diag(ones(3,1),-1);
% A = kron(eye(4),T) + diag(ones((4)^2  - 4,1),4) + diag(ones((4)^2 -4,1),-(4));
% u_1 = @(x,y) x.*(x-1);
% h = 1/(n-1);
% x = 0:h:1;
% y = 0:h:1;
% [X,Y] = meshgrid(x,y); 
% sol = u_1(X,Y);
% sol = reshape(sol.', [],1);
% b = A*sol;
% size(b)
% d2u_1 = @(x,y) -2.*ones(size(x));
% r = b;
% tol = 1e-10;
% k = 1;
% x = zeros(16,1);
% delta = r'*r;
% b_d = b'*b;
% p = r;
% 
% % while r'*r > tol
% for i = 1:4
%     s = twod_mult_ax(p, n,h);
%     a = delta/(p'*s);
%     x = x+a*p;
%     r = r-a*s;
%     if( r'*r < tol)
%         break;
%     end
%     delta_n = r'*r;
%     p = r + (delta_n/delta)*p;
%     delta = delta_n;
%     
% end