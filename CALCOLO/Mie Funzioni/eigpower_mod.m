function [lambdak,x,iter]=eigpower_mod(A,x0,tol,nmax)
%EIGPOWER    Approssima l'autovalore di modulo massimo
%            di una matrice.
%   LAMBDAK = EIGPOWER_MOD(A) calcola con il metodo delle
%   potenze l'autovalore di una matrice A di modulo
%   massimo a partire da un dato iniziale pari al
%   vettore unitario.
%   LAMBDAK = EIGPOWER_MOD(A,X0,TOL,NMAX) arresta il metodo
%   quando la differenza  fra due iterate consecutive
%   e' minore di TOL (il valore di default e' 1.E-06)
%   o quando il massimo numero di iterazioni NMAX (il
%   valore di default e' 100) e' stato raggiunto.
%   [LAMBDAK,X,ITER] = EIGPOWER_MOD(A,X0,TOL,NMAX)
%   restituisce anche l'autovettore unitario X tale
%   che A*X=LAMBDA*X ed il numero di iterazioni
%   effettuate per calcolare X.
%
%	LAMBDAK = vettore colonna contenente LAMBDA ad ogni iterazione

[n,m] = size(A);
if n ~= m
    error('Solo per matrici quadrate');
end

flag = '';
if nargin < 4
	nmax = 100;
end
if nargin < 3
	tol = 1.e-06;
end
if nargin < 2
	x0 = ones(n,1);
	flag = ', default vector';
end

fprintf("\nMetodo delle potenze (tol = %d, nnmax = %d%s)\n", tol, nmax, flag);

% iterazione zero fuori dal ciclo while
iter = 0;
y = x0/norm(x0); % y0
lambda = y'*A*y; % lambda0
err = tol*abs(lambda) + 1; % dobbiamo entrare nel ciclo
lambdak = lambda;
x = x0;
while err > tol*abs(lambda) && abs(lambda) ~= 0 && iter<nmax
   iter = iter + 1; % iter=1,2,3,...
   x = A*y;
   y= x/norm(x);
   lambdanew = y'*A*y;
   err = abs(lambdanew - lambda);
   lambda = lambdanew;
   lambdak = [lambdak; lambda];
end

if (err <= tol*abs(lambda))
     fprintf('Il metodo delle potenze converge in %d iterazioni  all''autovalore\n', iter);
	 lambda
else
     fprintf('Il metodo delle potenze non converge in %d iterazioni. \n', iter)
end

return
end