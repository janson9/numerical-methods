%% Bisection demo

f = @(x) x.^5 - 3*x + 1/2

% show more digits
format long g

a = -10;
b = 10;

assert (f(a)*f(b) < 0)

% exact answer (from where?)
exact = -1.35480588820312
clf;
xlabel('k')
ylabel('err')

tol = 1e-13

k = 0
while (b-a > tol)
  k = k + 1
  m = (a+b)/2;
  xk = m
  err = abs(xk - exact)
  if (f(m)*f(a) < 0)
    b = m;
  else
    a = m;
  end
  semilogy(k, err, 'ko')
  hold on;
  
end

