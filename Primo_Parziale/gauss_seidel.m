%%  GAUSS SEIDEL
%
%   [x,k,r,err] = gauss_seidel(A,b,x0,toll,nmax)
function [x,k,r,err] = gauss_seidel(A,b,x0,toll,nmax)
    n = length(b);
    if size(A,1) ~= n || size(A,2) ~= n ||  length(x0) ~= n
        error('Dimensioni sistema errate!');
    end
    P = tril(A);
    if prod(diag(P)) == 0
        error('Un elemento sulla diagonale della precondizionata è nullo');
    end
    x = x0;
    k = 0;
    r = b-A*x0;
    err = norm(r)/norm(b);
    while err > toll && k < nmax
        k = k+1;
        z = fwsub(P,r);
        x = x + z;
        r = b - A*x;
        err = norm(r)/norm(b);
    end
end