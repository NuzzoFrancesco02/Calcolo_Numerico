function [x,it] = ptofis_vettoriale(x0,phi,nmax,toll)
%
% [succ,it] = ptofis(x0,phi,nmax,toll,[a b]);
%
% --------Parametri di ingresso:
% x0      Punto di partenza
% phi     Funzione di punto fisso (definita inline o anonimous)
% nmax    Numero massimo di iterazioni
% toll    Tolleranza sul test d'arresto
% a b     Estremi dell'intervallo di esplorazione, SOLO PER OUTPUT GRAFICO
%
% --------Parametri di uscita:
% xvect   Vett. contenente tutte le iterate calcolate
%         (l'ultima componente e' la soluzione)
% it      Iterazioni effettuate
%

x = x0;
x = phi(x);
it =1;
err=toll+1;
while (err >= toll && it < nmax)
    xold = x;
    x = phi(x);
    err = norm(x-xold);
    it = it + 1;
end
end