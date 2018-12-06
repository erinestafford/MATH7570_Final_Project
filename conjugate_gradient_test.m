function x = conjugate_gradient_test(np, rhs,x0,x1,y0,y1)
num =np;
x = zeros(num+2,num+2); %phi
dx = 1/(np+1);
r = rhs;
p = r;
Ap = r;

rr = 0.0;
for i = 1:num
    for j = 1:num
        p(i,j) = r(i,j);
        rr = rr+ (r(i,j) * r(i,j));
    end
end

k = 0;
tol = 1.0e-8;
while (rr > tol) 
		for i = 2:num+1
            for j = 2:num+1
                Ap(i,j) = -4.0 * p(i,j)+ p(i-1,j)+p(i,j+1)+ p(i+1,j)+ p(i,j-1);
            end
        end

        pAp = 0.0;
		for i = 2:num+1
            for j = 2:num+1
                pAp = pAp+ p(i,j) * Ap(i,j);
            end
        end
		alpha = rr / pAp;

		rr1 = 0.0;
		for i = 1:num+2
            for j = 1:num+2
                %Apply BCs
            if(i ==1)
                x(i,j) = x0(j);
                r(i,j) = 0;
            elseif(i==np+2)
                x(i,j) = x1(j);
                r(i,j) = 0;
            elseif(j==1)
                x(i,j) = y0(i);
                r(i,j) = 0;
            elseif(j==np+2)
                x(i,j) = y1(i);
                r(i,j) = 0;
            else
                x(i,j) = x(i,j) + alpha * p(i,j);
                r(i,j) = r(i,j) - alpha * Ap(i,j);
                rr1 =rr1+ (r(i,j) * r(i,j));
            end
            end
        end
        
        beta = rr1 / rr;

		for i = 2:num+1
            for j = 2:num+1
                p(i,j) = r(i,j) + beta * p(i,j);
            end
        end

		rr = rr1;
		k = k+1;
	end

end