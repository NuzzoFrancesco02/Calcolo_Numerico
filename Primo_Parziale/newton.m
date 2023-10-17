%% [xvect, it] = newton(x0, toll, nmax, fun, dfun, mol)
% ATTENZIONE!!!     x(1) = x0
function [xvect, it] = newton(x0, toll, nmax, fun, dfun, mol)
if nargin == 5 || mol == 1
    mol = 1;
    fprintf('Lo zero ha molteplicità %d, convergerà quadraticamente \n', mol);
else
    fprintf('Lo zero ha molteplicità %d, convergerà linearmente \n', mol);
end
it = 0;
err = toll + 1;
x = x0;
xvect = x0;
while err > toll && it < nmax
    fc = fun(x);
    dfc = dfun(x);
    if (dfc == 0)
        error('La derivata prima è nulla alla %d° iterata \n', it);
    end
    x_vec = x;
    x = x - mol*fc/dfc;
    xvect = [xvect; x];
    err = abs(x-x_vec);
    it = it + 1;
end
    if it < nmax
        fprintf('Convergenza al passo %d \n', it);
    else
        fprintf('È stato raggiunto il numero massimo di iterazioni: %d \n', it);
    end
fprintf(' Radice calcolata : %-12.8f \n\n', xvect(end));
end