A = getA(4);
sol = [ 0
     0
     0
     0
     0
   1
   2
     0
     0
   2
   1
     0
     0
     0
     0
     0]; %(solution)
b = A*sol;

x = conjugate_gradient_solve(2, b,0,0,0,0);
x = x(:)
