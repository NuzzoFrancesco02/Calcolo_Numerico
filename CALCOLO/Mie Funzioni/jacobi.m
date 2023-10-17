%Metodo di jacobi

function[x,k] = jacobi(A,b,x0,tol,nmax)

n = length(b);
k = 0;

if size(A,1) ~= n || size(A,2) ~= n || length(x0) ~= n
    error('Le dimensioni non sono compatibili')
end

if prod(diag(A)) == 0
    error('Elemeni diagonali nulli')
end

D = diag(A);
x = x0;
r = b - A*x;
res_norm = norm(r)/norm(b);
while res_norm > tol && k < nmax
    k = k + 1;
    z = D \ r;
    x = x + z;
    r = b - A*x;
    res_norm = norm(r)/norm(b);
end

if res_norm < tol
	fprintf('Jacobi converge in %d iterazioni\n', k)
else
	fprintf('Jacobi non converge in %d iterazioni\n', nmax)
end
end