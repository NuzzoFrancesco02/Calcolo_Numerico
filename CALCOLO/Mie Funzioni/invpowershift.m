function [lambda,x,iter]=invpowershift(A,s,x0,tol,nmax)

%INVPOWERSHIFT  Approssima l'autovalore di una matrice
%               piu' vicino al numero (complesso) mu
%   LAMBDA = INVPOWERSHIFT(A,S,TOL,NMAX,X0) calcola con il metodo delle
%   potenze inverse con shift l'autovalore di una matrice A 
%   piu vicino a mu a partire da un dato iniziale pari
%   al vettore unitario.
%   LAMBDA = INVPOWERSHIFT(A,S,TOL,NMAX,X0) arresta il metodo
%   quando la differenza  fra due iterate consecutive
%   e' minore di TOL (il valore di default e' 1.E-06)
%   o quando il massimo numero di iterazioni NMAX (il
%   valore di default e' 100) e' stato raggiunto.
%   [LAMBDA,X,ITER] = INVPOWERSHIFT(A,TOL,NMAX,X0)
%   restituisce anche l'autovettore unitario X tale
%   che A*X=LAMBDA*X ed il numero di iterazioni
%   effettuate per calcolare X.

[n,m] = size(A);
if n ~= m, error('Solo per matrici quadrate'); end

flag = '';
if nargin < 5
	nmax = 100;
end
if nargin < 4
	tol = 1.e-06;
end
if nargin < 3
	x0 = ones(n,1);
	flag = ', default vector';
end

fprintf("\nMetodo delle potenze inverse con shift (tol = %d, nnmax = %d%s)\n", tol, nmax, flag);

M = A - s*eye(n);
% calcolo la fattorizzazione LU una volta per tutte
[L,U,P]=lu(M);
% iterazione zero fuori dal ciclo while
iter = 0;
y = x0/norm(x0); % y0
lambda = y'*A*y; % lambda0
err = tol*abs(lambda) + 1; % dobbiamo entrare nel ciclo
x = x0;
while (err>tol*abs(lambda)) && (abs(lambda)~=0) && (iter<nmax)
   iter = iter + 1; % iter=1,2,3,...
   % risolvo Mx^{(k)}=y^{(k-1)}
   z=fwsub(L,P*y);
   x=bksub(U,z);
   y= x/norm(x);
   lambdanew = y'*A*y;
   err = abs(lambdanew - lambda);
   lambda = lambdanew; 
end
if (iter < nmax)
%     fprintf('Il metodo delle potenze inverse con shift converge in %d iterazioni all''autovalore lambda = %d\n', iter, lambda);
     fprintf('Il metodo delle potenze inverse con shift converge in %d iterazioni all''autovalore \n', iter);
	lambda
else
     fprintf('Il metodo delle potenze inverse con shift non converge in %d iterazioni \n', iter);
end
return
end