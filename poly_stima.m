%% STIMA ERRORE POLINOMIO DI LAGRANGE
%       INPUT:
%   f : funzione
%   n : grado del polinomio
%   a, b : estremi dell'intervallo
%   str =>  'equi' : nodi equispaziati
%           'CGL' : nodi di Chebyshev–Gauss–Lobatto
%
function stim = poly_stima(f,n,a,b,str)
    der = f;
    x = a:0.001:b;
    for i = 1 : n + 1
         der = Jac(der);
    end

    if strcmp(str,'equi')
        h = (b-a)/n;
        x_nod = a:h:b;
        stim = max(abs(der(x)))*(((b-a)/n)^(n+1))/(4*(n+1));
    elseif strcmp(str,'CGL')
        i = 0:n;
        x_nod = -cos(pi.*i./n);
        syms y;
        omega = 1;
        for i = 1 : n + 1
            omega = omega*(y-x_nod(i));
        end
        omega = matlabFunction(simplify(expand(omega)));
        stim = max(abs(der(x)))*max(abs(omega(x)))/factorial(n+1);
    end