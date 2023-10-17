%% function I = integr_doppio(M1,M2,fun,str,a,b,c,d)
% str : 't' : trapezi
%       'pm_ret' : punti medi con rettangolo --> M1, M2 : num. lati NON num. rettangoli
%       'pm_tri' : punti medi con triangolo --> M1 : vettore RIGA
%                                                    contenente coordinate x
%                                               M2 : vettore RIGA
%                                                    contenente coordinate y
%       'GL_1' : formula di quadratura di Gauss-Legendre di ordine 1
function I = integr_doppio(M1,M2,fun,str,a,b,c,d)
    if nargin == 7
        M2 = M1;
    end
    if strcmp(str,'t')
        I = (b-a)*(d-c)*(fun(a,c)+fun(a,d)+fun(b,c)+fun(b,d))/4;
    elseif strcmp(str,'pm_ret') && (nargin == 8 || nargin == 7)
        h1 = (b-a)/M1;
        h2 = (d-c)/M2;
        x = a+h1/2:h1:b-h1/2;
        y = c+h2/2:h2:d-h2/2;
        S = 0;
        for k = 1 : M1
            for j = 1 : M1
                S = S + fun(x(k),y(j));
            end
        end
        I = (b-a)*(d-c)*S/(M1*M2);
    elseif strcmp(str,'pm_tri') && nargin == 4
        M = [M1' M2' ones(3,1)];
        A = 0.5*det(M);
        P = 0;
        for i = 1 : 3
            P = P + [M1(i) M2(i)];
        end
        P = P/3;
        I = A*fun(P(1),P(2));
    elseif strcmp(str,'GL_1')
        e = [-1/sqrt(3) 1/sqrt(3)];
        x = (a+b)/2 + ((b-a)/2).*e;
        y = (c+d)/2 + ((d-c)/2).*e;
        S = 0;
        for i = 1 : 2
            for j = 1 : 2
                S = S + fun(x(i),y(j));
            end
        end
        I = S*(b-a)*(d-c)/4;
    end
    