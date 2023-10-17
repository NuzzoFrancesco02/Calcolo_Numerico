%% x = gaussnewton(F,J,x0,toll,nmax)
function x = gaussnewton(F,J,x0,toll,nmax)
err = toll + 1;
k = 0;
x = x0;
while k < nmax && err > toll
    B = J(x);
    d = (B'*B)\(-B'*F(x));
    x = x + d;
    err = abs(max(F(x)));
    k = k + 1;
end
