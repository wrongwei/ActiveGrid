
c_eps = .5;
L = .1;
u = .14;
u_mean = 1.4;
nu = 1.5e-5;
eps = c_eps*u^3/L;
eta = nu^(3/4)*eps^(-1/4);
frq = u_mean/eta;
disp(eps);
disp(eta);
disp(frq);