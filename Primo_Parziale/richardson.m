function [x, k, r, alpha] = richardson(A, b, P, x0, toll, nmax, alpha)
% [x, k, r, alpha] = richardson(A, b, P, x0, toll, nmax, alpha)
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
% r: residuo
% alpha: angolo alpha(it-1) usato per calcolare x(it)
%
    n = length(b);
    if (size(A,1)~=n || size(A,2)~=n)
        error('Le dimensioni sono incompatibili');
    end

    x = x0;
    r = b - A*x;
    err = norm(r)/norm(b);
    k = 0;
    while err > toll && k < nmax
        k = k + 1;
        z = P \ r;
        if nargin == 6
            alpha = (z'*r)/(z'*A*z); % Dinamico
        end
        x = x + alpha*z;
        r = b - A*x; %r = r - aplha*A*z
        err = norm(r)/norm(b);
    end
end