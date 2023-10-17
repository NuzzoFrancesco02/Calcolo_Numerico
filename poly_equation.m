%% EQUAZIONE DEL POLINOMIO DI LAGRANGE
% [eq, phi, A, w] = poly_equation(n,a,b,str)
% Restituisce un vettore contente le phi(k)
function [eq, phi, A, w] = poly_equation(n,a,b,str)
    if strcmp(str,'equi')
        h = (b-a)/n;
        x_nod = a:h:b;
        syms x;
        phi = [];
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
        if nargin == 4
            P = simplify(expand(phi));
            S = 0;
            x = a:0.01:b;
            for i = 1 : n+1
                phi = matlabFunction(P(i));
                S = S + abs(phi(x));
            end
            A = max(S);
        end
        syms x;
        w = 1;
        for i = 0 : n
            w = (x-x_nod(i+1)).*w;
        end
    elseif strcmp(str,'CGL')
        k = 0:n;
        t = -cos(pi*k/n);
        x_nod = ((b-a)/2)*t + (a+b)/2;
        syms x;
        phi = [];
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
        if nargin == 4
            P = simplify(expand(phi));
            S = 0;
            x = a:0.01:b;
            for i = 1 : n+1
                phi = matlabFunction(P(i));
                S = S + abs(phi(x));
            end
            A = max(S);
        end
        syms x;
        w = 1;
        for i = 0 : n
            w = (x-x_nod(i+1)).*w;
        end
    end

    
    