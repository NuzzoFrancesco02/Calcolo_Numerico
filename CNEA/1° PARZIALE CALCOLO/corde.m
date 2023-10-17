% [x, it] = corde(f, a, b, x0, nmax, toll)
function [x, it] = corde(f, a, b, x0, nmax, toll)
    x = x0;
    q = (f(b)-f(a))/(b-a);
    it = 0;
    err = toll + 1;
    while err > toll && it < nmax
        x = x - (f(x)/q);
        err = abs(f(x));
        it = it + 1;
    end