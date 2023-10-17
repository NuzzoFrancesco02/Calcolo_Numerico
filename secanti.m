function [x,i,tolf]=secanti(x0,x1,f,tolx,nmax)
%SECANTI Esegue il metodo delle secanti, per la risoluzione di f(x)=0
%
%   [x,i,tolf]=SECANTI(x0,x1,f,tolx,nmax)
%
%   I parametri della funzione sono:
%       x0 -> il punto iniziale e prima approssimazione di x
%       x1 -> la seconda approssimazione della soluzione x
%       f -> funzione di cui valutare uno zero
%       tolx -> tolleranza per la radice
%       nmax -> limite superiore al numero di iterazioni
%
%   I valori di ritorno sono:
%       x -> la soluzione trovata
%       i -> il numero di iterazioni impiegate per ottenere la soluzione
%       tolf -> la tolleranza sulla funzione
%
%   See Also NEWTON, CORDE, STEFFENSEN
  i=0;
  fx0=feval(f,x0);
  err=abs(x1-x0);
  while (i<nmax & err>tolx)
      fx1=feval(f,x1);
      dfx1=(fx1-fx0)/(x1-x0);
      tolf=tolx*abs(dfx1);
      if abs(fx1)<=tolf
         break
      end
      x2=x1-(fx1/dfx1);
      err=abs(x2-x1);
      x0=x1;
      x1=x2;
      fx0=fx1;
      i=i+1;
  end
  x=x1;
