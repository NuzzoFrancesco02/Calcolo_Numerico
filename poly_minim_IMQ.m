%% MINIMI IMQ
% Ritorna un vettore contente i coefficienti per cui il problema ai
% quadrati risulta minimo, utilizza la risoluzione del sistema lineare con
% la matrice di Vandermonde
% ATTENZIONE!!! a è ordinato come: a0 a1 a2...
%       INPUT
%   m : esponente più alto
%   x,y : punti dati dal testo
function a = poly_minim_IMQ(m,x,y)
    if size(y,2) ~= 1
        y = y';
    end
    n = length(y)-1;
   
    for i = 1 : m+1
        V(:,i) = x.^(i-1);
    end

    A = V'*V;
    q = V'*y;
    a = A\q;
    %a = flip(a);