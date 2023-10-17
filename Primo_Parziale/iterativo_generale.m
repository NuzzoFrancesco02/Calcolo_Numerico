function x = iterativo_generale(B,g,x0,nmax)
x = x0;
k = 0;
while k < nmax
    x = (B*x) + g;
    k = k + 1;
end