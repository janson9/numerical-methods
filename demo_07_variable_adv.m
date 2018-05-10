% explicit solution of advection eqn  u_t + a(x) u_x, periodic BCs

% Grid and initial data:
h = .01;
x = (-1+h:h:1)';
a = 0.5*sin(2*pi*x) + 0.25;
k = .25*h/max(abs(a));
% initial data:
%v = abs(x-0.0) < 0.1;
v = cos(3*pi*x);

figure(1); clf; hold on;
plot(x,a, 'r-', 'linewidth', 2);
plt = plot(x,v,'k.-', 'linewidth',4);
legend('wavespeed', 'soln');
axis([-1 1 -1.01 1.01])
grid on

% Sparse matrices for finite differences
N = length(x);  e = ones(N,1);
Db = spdiags([-e  e], [-1 0], N, N);  % backward difference
Df = spdiags([-e  e], [0  1], N, N);  % forward difference
Db(1,end) = -1;
Df(end,1) = 1;
Db = 1/h * Db;
Df = 1/h * Df;

Tf = 5;
% adjust either final time or time-step to have integer steps
numsteps = ceil(Tf / k);
%k = Tf / numsteps
Tf = k*numsteps;

disp('paused'); pause

% Time-stepping:
for n=1:numsteps
  % Note "vectorized" upwinding.
  vnew = v + k*(-a.*(...
      (a > 0).*(Db*v) + (a <= 0).*(Df*v) ));
  % Note: should think about how to do this with loops...
  % if (a(j) > 0)
  %    % use backward diff
  % else
  %    % use forward diff
  % end
  v = vnew;
  set(plt,'ydata',v)
  drawnow
  %pause
end

