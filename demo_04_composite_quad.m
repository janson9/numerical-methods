%% Composite integration demo

f = @(x) sin(x);

a = 0;
b = pi;
n = 2*1;  % even

% (we generally don't know this!)
Iexact = 2;

h = (b-a)/n;

fa = f(a);
fb = f(b);

%% trap
xmid = a + h*(1:(n-1));
I_trap = h/2 * (fa + fb + 2*sum(f(xmid)))

%% simpson
xeven = a + h*(2:2:(n-1));
xodd = a + h*(1:2:(n-1));
f_odd = f(xodd);
f_even = f(xeven);
I_simp = h/3 * sum([fa fb 4*f(xodd) 2*f(xeven)])

err_trap = I_trap - Iexact
err_simp = I_simp - Iexact
