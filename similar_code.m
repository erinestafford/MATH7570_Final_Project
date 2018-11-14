
%% For finite difference approach:

% Iteration part...

% Initialize solution (U) given iniital average concentrations:
% (I don't normalize to exactly the set average, but on average
% it will be that average:-)


U = zeros(nx,nx,nf); %concentration of each lipid type
% first two lipid types are inner leaflet, second two are outer
for i=1:nf
    U(:,:,i) = 2*rand(nx,nx,1)*conc(i); %random noise added to initial concentration
end

t0id = tic;

for istep=1:nsteps
  t1id = tic;
  %% Fix ghost bondaries on U (solution is periodic)
  U(1:nkern,:,:) = U(nx0+(1:nkern),:,:);
  U(nx-nkern+(1:nkern),:,:) = U(nkern+(1:nkern),:,:);

  U(:,1:nkern,:) = U(:,nx0+(1:nkern),:);
  U(:,nx-nkern+(1:nkern),:) = U(:,nkern+(1:nkern),:);
  

  %% Transform solution with FFT
  Upadhat = fft2(U);
  % Apply convolution in FFT space:
  Zhat = 0*Upadhat;
  for i=1:nf
    for j=1:nf
        Zhat(:,:,j) = Zhat(:,:,j) + Upadhat(:,:,i) .* phipadhat(:,:,(i-1)*nf+j); %adding ni* f_mm/dn_i, which is constant 
    end
  end
  % Back transform to real space:
  Z2 = real(ifft2(Zhat));

  %% dealing with convolution
  
  if 0
      % For testing only -- sets convolution to zero in
      % padded region:
      Z2(:,1:nkern,:) = 0;
      Z2(1:nkern,:,:) = 0;
      Z2(:,(pad-nkern+1):pad,:) = 0;
      Z2((pad-nkern+1):pad,:,:) = 0;
  else
      % Apply same boundary conditions (periodic) for
      % Convolutions as solution, so we can take derivatives
      % (central differences):
      Z2(1:nkern,:,:) = Z2(nx0+(1:nkern),:,:);
      Z2(nx-nkern+(1:nkern),:,:) = Z2(nkern+(1:nkern),:,:);
      Z2(:,1:nkern,:) = Z2(:,nx0+(1:nkern),:);
      Z2(:,nx-nkern+(1:nkern),:) = Z2(:,nkern+(1:nkern),:);
  end
  t1 = toc(t1id);
 %% Updating lipid concentrations - need to understand this code better
 
  % Index set to simplify numerical differentiation:
  jj = 2:nx-1;

  % Compute graidients and Laplacians using central differences:
  d2u = (1/dx^2) * (U(jj-1,jj,:)+U(jj+1,jj,:)+U(jj,jj-1,:,:)+U(jj,jj+1,:)-4*U(jj,jj,:));
  Z2x = (1/(2*dx)) * (Z2(jj,jj+1,:) - Z2(jj,jj-1,:));
  Z2y = (1/(2*dx)) * (Z2(jj+1,jj,:) - Z2(jj-1,jj,:));

  Ux = (1/(2*dx)) * (U(jj,jj+1,:) - U(jj,jj-1,:));
  Uy = (1/(2*dx)) * (U(jj+1,jj,:) - U(jj-1,jj,:));

  d2Z2 = (1/dx^2) * (Z2(jj-1,jj,:,:)+Z2(jj+1,jj,:,:)+Z2(jj,jj-1,:,:,:)+Z2(jj,jj+1,:,:)-4*Z2(jj,jj,:,:));

  % Assemble time derivative of solution using computed gradients
  % and Laplacians:
  Ut = zeros(nx,nx,nf);
  Ut(jj,jj,:) = d2u - (Ux.*Z2x + Uy.*Z2y + U(jj,jj,:).*d2Z2);
  
  %% noisescale = 0*diffcoeff;
  for i=1:nf
    %% noise = randn(nx,nx);
    %% noise = sqrt(dt) * (noise(jj+1,jj+1) + noise(jj+1,jj-1) + noise(jj-1,jj+1) + noise(jj-1,jj-1) - 4*noise(jj,jj));

    Ut(:,:,i) = diffcoeff(i) * Ut(:,:,i); %% + noisescale(i)*noise;
  end
  
  
  % Time integration step using forward Euler:
  U = U + dt*Ut;
  

  if istep == 1
      % Print time of first convolution operation (second or later
      % is probably more representative of average...):
    disp(sprintf('FFT    convolution time = %.3f ms',t1*1e3));
  end


end

% Print how long a timestep takes on average:
t0 = toc(t0id);
disp(sprintf('Integration time = %.3f s (%.3f ms/dt)',t0,t0*1e3/nsteps));

