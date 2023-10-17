function [x, k] = richardson(A, b, P, x0, tol, nmax, alpha)
% [x, k] = richardson(A, b, P, x0, toll, nmax, alpha)
%
% Metodo di Richardson stazionario precondizionato
%               o dinamico precondizionato (gradiente precondizionato)
%
% Parametri di ingresso:
%
% A: matrice di sistema
% b: vettore termine noto
% P: precondizionatore
% x: guess iniziale
% tol:  tolleranza criterio d'arresto
% nmax: numero massimo di iterazioni ammesse
% alpha: parametro di accelerazione; se non assegnato si considera 
%        il metodo dinamico (gradiente precondizionato)
%
% Parametri di uscita:
%
% x: soluzione
% k: numero di iterazioni
%

n = length(b);
if ((size(A,1) ~= n) || (size(A,2) ~= n) || (length(x0) ~= n))
  error('Dimensioni incompatibili')
end

% E' possibile utilizzare una sola variabile x al posto di xn e xv viste
% nel laboratorio precedente.

x = x0;
k = 0;
r    = b - A * x;
res_normalizzato  = norm(r) / norm(b);

while ((res_normalizzato > tol) && (k < nmax))
     z = P \ r;
     if (nargin == 6)
         alpha = (z' * r) / (z' * A * z); % alpha dinamico
     end
     x = x + alpha * z;
     r = b - A * x; % equivalente a: r = r - alpha * A * z;
     res_normalizzato  = norm(r) / norm(b);
     k = k + 1;
end

if (res_normalizzato < tol)
     fprintf('Richardson converge in %d iterazioni \n', k);
else
     fprintf('Richardson non converge in %d iterazioni. \n', k)
end
