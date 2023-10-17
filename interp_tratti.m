%% INTERPOLAZIONE A TRATTI DI GRADO n
%
%   [y0, y, x] = interp_tratti(f,a,b,n,int,str,x0)
%
%       INPUT:
%   f : funzione
%   a, b : estremi intervallo
%   n : grado
%   int : numero di intervalli
%   str =>  'equi' : nodi equispaziati
%           'CGL' : nodi di Chebyshev–Gauss–Lobatto
%   x0 : se si vuole il valore del polinomio del punto
%       OUTPUT:
%   y : punti del polinomio valutato nell'intervallo tra a e b
%   x : intervallo in cui viente valutato il polinomio
%   y0 : valore del polinomio valutato in x0

function [y0, y, x] = interp_tratti(f,a,b,n,int,str,x0)
    if strcmp(str,'equi')
        h = (b-a)/int;
        x_nod = a:h:b;
    elseif strcmp(str,'CGL')
        i = 0:int;
        x_nod = -cos(pi.*i./int);
    end
    y = [];
    x = [];
    x(1) = NaN;
    y(1) = NaN;
    
    for i = 1 : int
        j = (x_nod(i+1)-x_nod(i))/n;
        nods = x_nod(i):j:x_nod(i+1);
        x_cycl = x_nod(i):1e-5:x_nod(i+1);
        P = polyfit(nods,f(nods),n);
        y_pol = polyval(P,x_cycl);
        if  i ~= 1 || i ~= int
            x(end) = [];
            y(end) = [];
        end
        y = [y y_pol];
        x = [x x_cycl];
        
        if x0 >= x_nod(i) && x0 <= x_nod(i+1) && nargin == 7
            y0 = polyval(P,x0);
        end
    end
    
    
