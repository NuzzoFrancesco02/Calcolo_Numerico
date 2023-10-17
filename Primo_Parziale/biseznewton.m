function [xvect,it] = biseznewton(a,b,nmax_b,nmax_n,toll,fun,dfun,mol)
if nargin == 7
    mol = 1;
end
toll_bi = (b- a)/(2^(nmax_b+1));
[xvect_bi, it_bi] = bisez(a,b,toll_bi,fun);
[xvect_new, it_new] = newton(xvect_bi(end),toll,nmax_n,fun,dfun,mol);
xvect = [xvect_bi; xvect_new(2:end)];
it = it_bi + it_new;
%{
if fun(a)*fun(b) > 0
    error('La funzione deve avere segno opposto nei due estremi \n');
end
xvect = [];
it_bi = -1;
err = toll + 1;
nmax = ceil(log2((b-a)/toll)-1);
fprintf('\nMassimo numero di iterazioni ammissibili: %d\n', nmax);
    while err > toll && it_bi < nmax_b
        it_bi = it_bi + 1;
        x = (b+a)/2;
        fc = fun(x);
        if fc == 0 
            err = 0;
        else
            err = abs(fc);
        end
        xvect = [xvect;x];
        x0 = x; 
        if (fun(a)*fc) > 0
            a = x;
        else
            b = x;
        end
    end
if nargin == 7
    mol = 1;
end
it_new = 0;
err = toll + 1;
x = x0;
    while err > toll && it_new < nmax_n
        fc = fun(x);
        dfc = dfun(x);
        if (dfc == 0)
            error('La derivata prima è nulla alla %d° iterata \n', it_new);
        end
        x_vec = x;
        x = x - mol*fc/dfc;
        xvect = [xvect; x];
        err = abs(x-x_vec);
        it_new = it_new + 1;
    end
it = it_bi + it_new;
fprintf(' Radice calcolata : %-12.8f \n\n', xvect(end));
end
%}