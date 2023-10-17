%Metodo di Richardson

function[x,k] = richardson(A,b,P,x0,tol,nmax,alpha)

n = length(b);

if size(A,1) ~= n || size(A,2) ~= n || length(x0) ~= n
    error('Le dimensioni non sono compatibili')
end

if prod(diag(A)) == 0
    error('Elemeni diagonali nulli')
end

x = x0;
r = b - A*x;
res_norm = norm(r)/norm(b);
%B = eye(n) - alpha * (P^-1) * A;
k = 0;

while res_norm > tol && k < nmax
	k = k + 1;
	z = P \ r;
	if nargin == 6
		alpha = (z'*r)/(z'*A*z);
	end

	x = x + alpha*z;
%	r = r - alpha*A*z;
	r = b - A*x;
	res_norm = norm(r)/norm(b);
end

if res_norm < tol
	fprintf('Richardson converge in %d iterazioni\n', k)
else
	fprintf('Richardson non converge in %d iterazioni\n', nmax)
end
end