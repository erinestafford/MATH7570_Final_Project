function Ax = twod_mult_ax(x, n,h)
Ax = zeros(n, n);
x = reshape(x.', [(n), (n)])';
  for i = 1:n %x
        for j = 1:n %y
            Ax(i,j) = 0;
            % Apply BCs
            %corners
            if (i ==1 && j ==1) %x(1,1)
                Ax(i,j) = -4*x(i,j) + x(i,j+1) + x(i+1,j);
            elseif (i == 1 && j==n )% (1,n)
                Ax(i,j) = -4*x(i,j) + x(i+1,j) + x(i,j-1);
            elseif ( i == n && j==1 )% (n,1)
                Ax(i,j) = -4*x(i,j) + x(i-1,j) + x(i,j+1);
            elseif ( i == n && j==n )% (n,n)
                Ax(i,j) = -4*x(i,j) + x(i-1,j) + x(i,j-1);
            %bounds not corners
            elseif (i~=1 && i~=n && j==1) %column 1
                Ax(i,j) = -4*x(i,j) + x(i-1, j) + x(i,j+1) + x(i+1,j);
            elseif (i == 1 && j ~=n && j ~= 1)%row 1
                Ax(i,j) = -4*x(i,j) + x(i,j-1) + x(i,j+1) + x(i+1,j);
            elseif (i ~=1 && i ~= n && j==n) %column n
                Ax(i,j) = -4*x(i,j) + x(i-1, j) + x(i,j-1) + x(i+1,j);
            elseif (i == n && j~=n&&j ~=1) %row n
                Ax(i,j) = -4*x(i,j) + x(i,j-1) +x(i,j+1) + x(i-1,j);
            else
                Ax(i,j) = x(i-1,j)+x(i+1,j) + x(i,j-1) + x(i,j+1)-4*x(i,j);
            end
        end
  end
   Ax = Ax;
   Ax = reshape(Ax.', [],1);
end