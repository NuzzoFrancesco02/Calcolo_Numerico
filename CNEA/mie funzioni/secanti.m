function [x, it] = secanti(f, x0, q0, nmax, toll)
    x = x0;
    q = q0;
    it = 0;
    err = toll + 1;
    while err > toll && it < nmax
        xv = x;
        x = x - f(x)/q;
        q = (f(x)-f(xv))/(x-xv);
        err = abs(f(x));
        it = it + 1;
    end