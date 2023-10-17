%Metodo di Gauss-Seidel

function[x,k] = gauss_seidel(A,b,x0,tol,nmax)

n = length(b);
k = 0;

if size(A,1) ~= n || size(A,2) ~= n || length(x0) ~= n
    error('Le dimensioni non sono compatibili')
end

if prod(diag(A)) == 0
    error('Elemeni diagonali nulli')
end

T = tril(A);
x = x0;
r = b - A*x;
res_norm = norm(r)/norm(b);
while res_norm > tol && k < nmax
    k = k + 1;
    z = T \ r;
    x = x + z;
    r = b - A*x;
    res_norm = norm(r)/norm(b);
end

if res_norm < tol
	fprintf('Gauss-Seidel converge in %d iterazioni\n', k)
else
	fprintf('Gauss-Seidel non converge in %d iterazioni\n', nmax)
end
end