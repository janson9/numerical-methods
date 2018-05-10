%% heat equation approximate solution
% u_t = u_xx
% on x in [0, 1] with initial condition
% u(x, 0) = sin(2*pi*x) and periodic boundary
% conditions

N = 20;
h = 1/N;
x = [0:h:(1-h)]';   % later: why 1-h?

e = ones(N, 1);

% make a "sparse matrix":
L = spdiags([e -2*e e], [-1 0 1], N, N);
L(1, N) = 1;
L(N, 1) = 1;

%full(L) to see the matrix

% initial condition
u0 = sin(2*pi*x);

clf;
lw = 'linewidth';
plot(x, u0, 'rx--', lw, 3);
hold on;
H = plot(x, u0, 'bx-', lw, 3);
xlabel('x'); ylabel('u');

k = h^2/4
%k = h^2/1.5   % did not work

u = u0;

Tf = 0.125;
numsteps = ceil(Tf/k);

for n=1:numsteps
  % formula from the board, rearranged
  unew = u + k/h^2*(L*u);
  u = unew;
  set(H, 'ydata', u);
  drawnow;
end
