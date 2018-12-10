function sol = compute_gridpoints_fns(fn,x,y)
for i = 1:length(x)
    for j = 1:length(y)
        sol(i,j) = fn(x(i),y(j));
    end
end
end