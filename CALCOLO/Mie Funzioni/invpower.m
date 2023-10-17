function [lambda,x,iter]=invpower(A,x0,tol,nmax)
%INVPOWER    Approssima l'autovalore di modulo minimo
%            di una matrice.
%   [LAMBDA,X,ITER] = INVPOWER(A, X0, TOL, NMAX) 
%   calcola con il metodo delle potenze inverse l'autovalore 
%   di una matrice A di modulo minimo a partire da un dato 
%   iniziale pari al vettore unitario, arrestando il metodo
%   quando la differenza  fra due iterate consecutive
%   e' minore di TOL (il valore di default e' 1.E-06)
%   o quando il massimo numero di iterazioni NMAX (il
%   valore di default e' 100) e' stato raggiunto.
%   La funzione restituisce anche l'autovettore unitario 
%   X tale che A*X=LAMBDA*X ed il numero di iterazioni
%   effettuate per calcolare X.

[n,m] = size(A);
if n ~= m, error('Solo per matrici quadrate'); end

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

fprintf("\nMetodo delle potenze inverse (tol = %d, nnmax = %d%s)\n", tol, nmax, flag);

% calcolo la fattorizzazione LU una volta per tutte
[L,U,P]=lu(A);

% iterazione zero fuori dal ciclo while
iter = 0;
y = x0/norm(x0); % y0
lambda = y'*A*y; % lambda0
err = tol*abs(lambda) + 1; % dobbiamo entrare nel ciclo
x = x0;
while err > tol*abs(lambda) && abs(lambda) ~= 0 && iter<nmax
   iter = iter + 1; % iter=1,2,3,...
   % risolvo Ax^{(k)}=y^{(k-1)}
   z=fwsub(L,P*y);
   x=bksub(U,z);
   y= x/norm(x);
   lambdanew = y'*A*y;
   err = abs(lambdanew - lambda);
   lambda = lambdanew; 
end

if (err <= tol*abs(lambda))
     fprintf('Il metodo delle potenze inverse converge in %d iterazioni all''autovalore\n', iter);
	 lambda
else
     fprintf('Il metodo delle potenze inverse non converge in %d iterazioni \n', iter);
end

return
end