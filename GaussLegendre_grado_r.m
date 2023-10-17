function r = GaussLegendre_grado_r(a,b,h)
n = (b-a)/h-2
if mod(n,2)==0
    r = n + 1;
else 
    r = n;
end