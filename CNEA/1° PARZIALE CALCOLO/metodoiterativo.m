function [x,k] = metodoiterativo(A,b,x0,toll,nmax,w)
n = length(b);
    if size(A,1) ~= n || size(A,2) ~= n ||  length(x0) ~= n
        error('Dimensioni sistema errate!');
    end
D = diag(diag(A));
E = -tril(A,-1);
F = -triu(A,1);
B = (1-w)*eye(n) + w*((D-E)\F);
g = w*((D-E)\b);
x = x0;
k = 0;
r = b-A*x0;
err = norm(r)/norm(b);
    while err > toll && k < nmax
        k = k + 1;
        x = B * x + g;
        r = b - A*x;
        err = norm(r)/norm(b);
    end
