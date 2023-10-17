function [lambda,x,iter]=invpower(A,tol,nmax,x0)
%INVPOWER    Approssima l'autovalore di modulo minimo
%            di una matrice.
%   [LAMBDA,X,ITER] = INVPOWER(A, TOL, NMAX, X0) calcola con 
%   il metodo delle potenze inverse l'autovalore di una 
%   matrice A di modulo minimo a partire da un dato iniziale pari
%   al vettore unitario, arrestando il metodo
%   quando la differenza  fra due iterate consecutive
%   e' minore di TOL (il valore di default e' 1.E-06)
%   o quando il massimo numero di iterazioni NMAX (il
%   valore di default e' 100) e' stato raggiunto.
%   La funzione restituisce anche l'autovettore unitario 
%   X tale che A*X=LAMBDA*X ed il numero di iterazioni
%   effettuate per calcolare X.

[n,m] = size(A);
if n ~= m, error('Solo per matrici quadrate'); end
if nargin == 1
   tol = 1.e-06;   x0 = ones(n,1);   nmax = 100;
end

% calcolo la fattorizzazione LU una volta per tutte
[L,U,P]=lu(A);

% iterazione zero fuori dal ciclo while
iter = 0;
y = x0/norm(x0); % y0
lambda = y'*A*y; % lambda0
err = tol*abs(lambda) + 1; % dobbiamo entrare nel ciclo

while (err>tol*abs(lambda)) && (abs(lambda)~=0) && (iter<nmax)
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
     fprintf('Il metodo delle potenze inverse converge in %d iterazioni \n', iter);
else
     fprintf('Il metodo delle potenze inverse non converge in %d iterazioni. \n', iter)
end

return
