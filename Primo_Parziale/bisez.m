%% [xvect,it] = bisez(a,b,toll,fun)
% Parte da x(0)
%
% k_min > log2((b-a)/tol)-1
function [xvect,it] = bisez(a,b,toll,fun)

if fun(a)*fun(b) > 0
    error('La funzione deve avere segno opposto nei due estremi \n');
end

it = -1;
xvect = [];
err = toll + 1;
nmax = ceil(log2((b-a)/toll)-1);
fprintf('\nMassimo numero di iterazioni ammissibili: %d\n', nmax);

    while err > toll && it < nmax
        it = it + 1;
        x = (b+a)/2;
        fc = fun(x);
        if fc == 0 
            err = 0;
        else
            err = abs(fc);
        end
        
        xvect = [ xvect; x]; 

        if (fun(a)*fc) > 0
            a = x;
        else
            b = x;
        end
    end
    if (it == nmax)
        fprintf('Max numero di iterazioni raggiunto! Errore %6.4e \n', err);
    else 
        fprintf('x_%d soddisfa la tolleranza sul residuo \n', it);
    end
    
    fprintf('Radice calcolata: %-12.8f \n', xvect(end));
end