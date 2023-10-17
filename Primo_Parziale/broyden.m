function [x, B, k] = broyden(F,B0,x0,tol,nmax)
x = x0;
B = B0;
err = tol+1;
k = 0;
while err > tol && k < nmax
    d = B\(-F(x));
    xv = x;
    x = x + d;
    B = B + ((F(x)-F(xv)-B*d)*d')./(d'*d);
    k = k + 1;
    err = norm(F(x));
end