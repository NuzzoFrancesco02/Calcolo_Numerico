function [x,it]=jacobi(A,b,x0,toll,nmax)

% Metodo di Jacobi
%
% A: matrice del sistema
% b: termine noto
% x0: vettore iniziale
% toll: tolleranza sul residuo normalizzato
% nmax: massimo numero di iterazioni
%
% x: soluzione ottenuta
% it: numero di iterazioni effettuate

n = length(b);
it = 0;

if ((size(A,1) ~= n) || (size(A,2) ~= n) || (length(x0) ~= n))
  error('Dimensioni incompatibili')
end

if (prod(diag(A)) == 0)
  error('Errore: elementi diagonali nulli')
end

D_inv = diag(1./diag(A));
x = x0;
r = b - A*x;
err = norm(r) / norm(b);

while (err > toll && it < nmax)
    it = it + 1;
    z = D_inv*r;
    x = x + z;
    r = b - A*x;
    err = norm(r)/norm(b);
end



