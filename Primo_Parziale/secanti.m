%% [xvect, it] = secanti(f, x0,x1, nmax, toll)
% Vettore parte da 0, per avere x(2)-->x(2+1)
function [xvect, it] = secanti(f, x0,x1, nmax, toll)
    xvect(1) = x0;
    xvect(2) = x1;
    it = 1;
    err = toll + 1;
    while err > toll && it < nmax
        it = it + 1;
        q = (f(xvect(it))-f(xvect(it-1)))/(xvect(it)-xvect(it-1));
        xvect(it+1) = xvect(it) - f(xvect(it))/q;
        err = abs(f(xvect(end)));
    end