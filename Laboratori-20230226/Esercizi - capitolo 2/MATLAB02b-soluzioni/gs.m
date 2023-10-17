function [x,k]=gs(A,b,x0,toll,nmax)

% Metodo di Gauss-Seidel
%
% A: matrice del sistema
% b: termine noto
% x0: vettore iniziale
% toll: tolleranza sul residuo normalizzato
% nmax: massimo numero di iterazioni
%
% xn: soluzione ottenuta
% k: numero di iterazioni effettuate

n = length(b);
xn = zeros( n, 1 );
k = 0;

if (( size(A,1)~=n) || (size(A,2)~=n) || (length(x0) ~= n) )
  error('dimensioni incompatibili')
end

if (prod(diag(A)) == 0)
    error('errore: elementi diagonali nulli')
end

T = tril(A);
xv = x0;
r = b - A * x0;
err = norm(r) / norm(b);

while ( err > toll && k < nmax )
  k = k + 1;
  z = fwsub(T,r);
  xn = xv + z;
  r = b - A*xn;  
  err = norm(r) / norm(b);
  xv = xn;
end

x = xn;

