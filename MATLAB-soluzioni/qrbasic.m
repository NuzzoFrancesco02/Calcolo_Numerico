function D=qrbasic(A,tol,nmax)
%  QRBASIC calcola gli autovalori di una matrice
%  D=QRBASIC(A,TOL,NMAX) calcola con il metodo delle
%  iterazioni QR tutti gli autovalori della matrice A
%  a meno di una tolleranza TOL ed in un numero massimo
%  di iterazioni NMAX. La convergenza di questo metodo
%  non e' in generale garantita.

[n,m]=size(A);
if n ~= m
  error('La matrice deve essere quadrata')
end
T = A; 
niter = 0; 
test = max(max(abs(tril(T,-1))));
while niter <= nmax && test >= tol
  [Q,R]=qr(T);    
  T = R*Q;
  niter = niter + 1;
  test = max(max(abs(tril(T,-1))));
end
if niter > nmax
 fprintf(['Il metodo non converge nel massimo',...
             ' numero di iterazioni permesso']);
else
 fprintf('Il metodo converge in %d iterazioni\n',niter)
end
D = diag(T);
return
