function K = cond_spect_sp(A,P,x0,toll,nmax)
x = x0;
y_max = x./norm(x);
y_min = y_max;
K = 1;
it = 0;
err = toll + 1;
Ra = chol(A);
Rp = chol(P);
while err > toll && it < nmax
    fp = fwsub(Rp',A*y_max);
    x_max = bksub(Rp,fp);
    y_max = x_max./norm(x_max);

    fa = fwsub(Ra',P*y_min);
    x_min = bksub(Ra,fa);
    y_min = x_min./norm(x_min);

    Kv = K;
    K = ((y_max'*A*y_max)/(y_max'*P*y_max))*((y_min'*P*y_min)/(y_min'*A*y_min));
    err = abs(Kv-K);
    it = it + 1;
end