function [xk,k] = richardson_it(A,b,P,x0,toll,nmax,alpha)

% Metodo del gradiente
%
% Parametri di ingresso:
%
% A: matrice di sistema
% b: vettore termine noto
% P: precondizionatore del metodo. Si usa P=I nel caso non precondizionato.
% x0: guess iniziale
% toll:  tolleranza criterio d'arresto
% nmax: numero massimo di iterazioni ammesse
% alpha: parametro di accelerazione; se non assegnato si considera 
%        il metodo dinamico (gradiente precondizionato)
%
% Parametri di uscita:
%
% xk: matrice le cui colonne sono i vettori della soluzione ad ogni iterazione
% k: numero di iterazioni

n = length(b);
if ((size(A,1) ~= n) || (size(A,2) ~= n) || (length(x0) ~= n))
  error('Dimensioni incompatibili')
end

x = x0;
k = 0;
r    = b - A * x;
res  = norm(r) / norm(b);
xk = [x0];

while ((res > toll) && (k < nmax))
     z = P \ r;
     if (nargin == 6)
         alpha = (z' * r) / (z' * A * z); % alpha dinamico
     end
     x = x + alpha * z;
     r = b - A * x;
     res  = norm(r) / norm(b);
     k = k + 1;
     xk=[xk x];
end

if (res < toll)
     fprintf('Richardson converge in %d iterazioni \n', k);
else
     fprintf('Richardson non converge in %d iterazioni. \n', k)
end

end