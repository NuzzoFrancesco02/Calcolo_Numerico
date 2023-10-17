%% POLINOMIO DI LAGRANGE
% INPUT
%   n : grado del polinomio
%   f : funzione esatta
%   a,b : estremi funzione
%   x : intervallo discreto
%   str =>  'equi' : nodi equispaziati
%           'CGL' : nodi di Chebyshev–Gauss–Lobatto
function [y_pol, err] = poly_lagr(f,n,a,b,x,str)
if strcmp(str,'equi')
    h = (b-a)/n;
    x_nod = a:h:b;
    y_pol = polyval(polyfit(x_nod,f(x_nod),n),x);
    err = abs(f(x)-y_pol);
elseif strcmp(str,'CGL')
    k = 0:n;
    t = -cos(pi*k/n);
    x_nod = ((b-a)/2)*t + (a+b)/2;
    y_pol = polyval(polyfit(x_nod,f(x_nod),n),x);
    err = abs(f(x)-y_pol);
end