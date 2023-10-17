%% SCARTO QUADRATICO
%   scarto = poly_scarto_IMQ(x,y,n)
%       INPUT:
%   x, y : punti
%   n : grado polinomio
%
function scarto = poly_scarto_IMQ(x,y,n)
    p2 = polyval(polyfit(x,y,n),x);
    scarto = sum((p2-y).^2);