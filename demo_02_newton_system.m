%% Newton demo
% system example

f = @(x) [x(1).*x(2) + x(2).^2 - 2;
          x(1).^3*x(2) - 3*x(1) - 1];

% Jacobian

% N^2 calculus 100 students ;-)
J = @(x) [x(2) ,   x(1) + 2*x(2);
          3*x(1).^2.*x(2) - 3 ,  x(1).^3];

x0 = [0; 1];

xk = x0;

% why not?
%tol = 1e-25;

tol = 1e-10

% big dummy value
delta = [inf; inf];

while (norm(delta) > tol)
  % backslash to solve J*delta = -f
  % (70 years of work to help us)
  delta = J(xk) \ -f(xk);
  xk = xk + delta;
  disp(xk')
end

% check!
f(xk)
