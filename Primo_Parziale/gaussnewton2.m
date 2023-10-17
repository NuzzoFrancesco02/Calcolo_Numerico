function x = gaussnewton2(F,J,x0,toll,nmax)
x = x0;
err = toll + 1;
it = 0;
while err > toll && it < nmax
    B = J(x);
    d = (B'*B)\(-B'*F(x));
    x = x + d;
    it = it + 1;
    err = norm(F(x));
end