function [D,Q,R] = qrbasic(A,tol,nmax)
%  QRBASIC calcola gli autovalori di una matrice
%  [D,Q,R]=QRBASIC(A,TOL,NMAX) calcola con il metodo delle
%  iterazioni QR tutti gli autovalori della matrice A
%  a meno di una tolleranza TOL ed in un numero massimo
%  di iterazioni NMAX. La convergenza di questo metodo
%  non e' in generale garantita.

[n,m]=size(A);
if n ~= m
	error('La matrice deve essere quadrata')
end

if nargin < 3
	nmax = 100;
end
if nargin < 2
	tol = 1.e-06;
end
fprintf("\nMetodo delle potenze inverse (tol = %d, nnmax = %d)", tol, nmax);

T = A; 
iter = 0; 
test = max(max(abs(tril(T,-1))));
while iter < nmax && test > tol
	[Q,R] = qr(T);    
	T = R*Q;
	iter = iter + 1;
	test = max(max(abs(tril(T,-1))));
end

if test > tol
	fprintf('\nIl metodo delle iterazioni QR non converge in %d iterazioni\n', nmax);
else
	fprintf('\nIl metodo delle iterazioni QR converge in %d iterazioni\n', iter)
end

D = diag(T);

return
end