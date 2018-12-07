function Ax = mult_ax(x, np)
    T = diag(-4*ones(np+2,1)) + diag(ones(np+1,1),1)+ diag(ones(np+1,1),-1);
    A = kron(eye(np+2),T) + diag(ones((np+2)^2  - (np+3),1),np+3) + diag(ones((np+2)^2 -(np+3),1),-(np+3));
    
    Ax = A*x;
    
end