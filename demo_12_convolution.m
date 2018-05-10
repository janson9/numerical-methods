% Gray--Scott reaction-diffusion using spectral convolution-based
% time-stepping

% Gray--Scott parameters and reaction functions
Fp = 0.05;  kp = 0.065;  Du = .16/(14^2);  Dv = Du/3;
%Fp = 0.054;  kp = 0.063;  Du = .16/(14^2);  Dv = Du/3;
f = @(u,v) (-u.*v.*v + Fp*(1-u));
g = @(u,v) ( u.*v.*v - (Fp+kp)*v);

% Computational domain
wx = 4;      wy = 4;
Nx = 128;     Ny = 128;
dx = wx/Nx;  dy = wy/Ny;
x = -wx/2:dx:wx/2-dx;  y = -wy/2:dy:wy/2-dy;
[xx,yy] = meshgrid(x,y);

dt = 2;  Tf = 10000;  % time-step size and time final time

% Analytically derived heat kernels
k = fftshift( (2*pi/wx)*[-Nx/2:Nx/2-1] );
l = fftshift( (2*pi/wy)*[-Ny/2:Ny/2-1] );
[kx,ky] = meshgrid(k,l);  rr = kx.^2 + ky.^2;
kernel_u = exp(-rr*Du*dt);  kernel_v = exp(-rr*Dv*dt);

% IC: a perturbation around 0 and 1
pert = exp(-(2*(xx-0.5).^2 + 1.5*(yy+0.25).^2)) + ...
       exp(-(2*(xx+0.5).^2 + 1.5*(yy-0.40).^2));
u = 1 - pert;
v = 0 + 0.5*pert;

disp('Starting time steps');
for i = 1:(Tf/dt)
  t = dt*i;
  % Add reactions then convolute with heat kernel
  u = real(ifftn(kernel_u .* fftn( u + dt*f(u,v) )));
  v = real(ifftn(kernel_v .* fftn( v + dt*g(u,v) )));

  if (mod(i,50)==0)
    disp([i t])
    figure(1); clf;
    pcolor(x,y,u)
    shading flat
    title(sprintf('t=%g (step num %d)', t, i));
    colorbar
    axis equal
    figure(2); clf;
    pcolor(x,y,v)
    shading flat
    colorbar
    axis equal
    pause(0);
    drawnow;
  end
end
