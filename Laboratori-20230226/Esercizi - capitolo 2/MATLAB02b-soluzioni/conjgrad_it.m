function [xk,it] = conjgrad_it(A,b,x0,nmax,toll)
x = x0;
r = b - A*x0;
p = r;
xk=x0;

res_norm = norm(r) / norm(b);
it = 0;
while it < nmax && res_norm > toll
    it = it + 1;
    
    alpha = (p' * r) / (p' * A * p);
    x = x + alpha * p;
    r = r - alpha * A * p;
    beta = (p' * A * r) / (p' * A * p);
    p = r - beta * p;
    
    res_norm = norm(r) / norm(b);
    xk = [xk, x];
end