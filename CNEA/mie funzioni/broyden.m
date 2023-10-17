function [x, B] = broyden(F,B0,x0,toll,nmax)
k = 0;
x = x0;
n = size(x0);
B = B0;
err = toll + 1;
while err > toll && k < nmax
    d = B \ (-F(x));
    x_old = x;
    x = x + d;
    B = B + (F(x)-F(x_old)-B*d)*(d')/dot(d,d);
    k = k + 1;
end
