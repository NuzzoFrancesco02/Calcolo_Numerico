%% EQUAZIONE DEL POLINOMIO DI LAGRANGE
% [phi, eq] = poly_equation(x)
% Restituisce un vettore contente le phi(k)
function [eq, phi] = poly_equation(x_nod)
    syms x;
    phi = []
    n = length(x_nod)-1;
    for k = 0:n
        phik = 1;
        for i = [0:k-1, k+1:n]
            phik = phik*((x-x_nod(1+i))./(x_nod(1+k)-x_nod(1+i)));
        end
        simplify(phik);
        phi = [phi phik];
        if k > 0
            eq = phi(end).*phik;
        end
    end
    eq = (expand(eq));
    
    