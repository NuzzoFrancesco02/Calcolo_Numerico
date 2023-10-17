function [x,k,xk] = conjgrad_it(A,b,x0,tol,nmax)

n = length(b);

if size(A, 1) ~= n || size(A, 2) ~= n || length(x0) ~= n
	error('Le dimensioni non sono compatbili')
end


if prod(diag(A)) == 0
	error('Elementi diagonali nulli')
end

if eig(A) <= 0
	error('Matrice non dp')
end

k = 0;
x = x0;
xk = x;
r = b - A*x;
p = r;
res_norm = norm(r) / norm(b);

while res_norm > tol && k < nmax
	alpha = (p'*r) / (p'*A*p);
	x = x + alpha*p;
	r = r - alpha*A*p;
	beta = (p'*A*r) / (p'*A*p);
	p = r - beta*p;
	res_norm = norm(r) / norm(b);
	k = k + 1;
	xk = [xk x];
end

if res_norm < tol
	fprintf('Il metodo del gradiente coniugato converge in %d iterazioni\n', k)
else
	fprintf('Il metodo del gradiente coniugato non converge in %d iterazioni\n', nmax)
end
end