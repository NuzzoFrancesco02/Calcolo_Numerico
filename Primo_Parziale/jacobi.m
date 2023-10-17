% [x,k] = jacobi(A,b,x0,toll,nmax)
%
%x : ultimo elemento
%k : ultima di iterazione
%A : matrice
%b : sol
%x0 : prima soluzione
%toll: tolleranza
%nmax: massima iterazione

function [x,k] = jacobi(A,b,x0,toll,nmax)
    n = length(b);
    if (size(A,1)~=n || size(A,2)~=n || length(x0)~=n)
        error('Le dimensioni non sono compatibili!')
    end
    if prod(diag(A))==0
        error('Un elemento diagonale nullo')
    end
    
    D_inv = diag(1./diag(A));
    x = x0;
    r = b-A*x;
    err = norm(r)/norm(b);
    k = 0;
    while( err > toll && k <= nmax)
        k = k+1;
        z = D_inv * r;
        x = x + z;
        r = b - A*x;
        err = norm(r)/norm(b);
    end
end

