function [xk,x, k, r] = richardson_it(A, b, P, x0, toll, nmax, alpha)
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
% xk : matrice per colonne tutte le iterata di x
% x: soluzione
% k: numero di iterazioni
%
    n = length(b);
    if (size(A,1)~=n || size(A,2)~=n)
        error('Le dimensioni sono incompatibili');
    end

    x = x0;
    r = b - A*x;
    err = norm(r)/norm(b);
    k = 0;
    xk = [x0];
    while err > toll && k < nmax
        k = k + 1;
        z = P \ r;
        if nargin == 6
            alpha = (z'*r)/(z'*A*z); % Dinamico
        end
        x = x + alpha*z;
        r = r - alpha*A*z;
        err = norm(r)/norm(b);
        xk = [xk x];
    end
end