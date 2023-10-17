function [t, u] = Runge_Kutta_gen(A,b,c,f,t_max,h,y0)
    s = size(A,1);
    if size(A,2) ~= s 
        error('Matrice A sbagliata');
    elseif length(c) ~= s
        error('Vettore c sbagliato');
    elseif length(b) ~= s
        error('Vettore b sbagliato');
    end
    if nargin == 7
        t = 0:h:t_max;
        N = t_max/h;
        u = [y0];
    if A ~= tril(A)
        error('A non triangolare inferiore');
    end
    for n = 1 : N
        K = [];
        for i = 1 : s
            S = 0;
            for j = 1 : s
                if j<i
                    S = S + A(i,j)*K(j);
                end
            end
            K = [K f(t(n)+c(i)*h,u(n)+h*S)];
        end
        u = [u u(end)+h*sum(b.*K)];
    end
    end
    if nargin >= 3
        if nargin ~= 7
            u = NaN;
        end
        syms h y0 l;
        assume(l<0);
        esp = 0:s;
        R = sum(((l*h).^esp)./(factorial(esp)));
        max = solve(abs(R)-1,h,"Real",true);
        max = double(max*l);
        fprintf('\n %.0f < h < %.4f/|l|\n\n',abs(max));
    end
    
    