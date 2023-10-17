%%
% N = numero intervalli. (Se semplice N = 1)
%
% n = num. nodi - 1
% str = 'equi','CGL'
function I = gausslegendre_comp(a,b,fun,N,n,str)
    I = 0;
    if strcmp(str,'equi')
    h = (b-a)/N;
    x = a:h:b;
        for i = 2 : N+1
            [xi, wi] = zplege(n, x(i-1), x(i));
            I = I + sum(wi.*fun(xi));
        end
    elseif strcmp(str,'CGL')
        i = 0:n;
        x = -cos(pi.*i./n);
        for i = 2 : N
            [xi, wi] = zplege(n, x(i-1), x(i));
            I = I + sum(wi.*fun(xi));
        end
    end