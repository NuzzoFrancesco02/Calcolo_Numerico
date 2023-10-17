%% Direzione ottimale al passo k-esimo del Gradiente Coniugato
% 
% [direct, teta] = grad_conj_direction(A,b,x0,max)
%   
%   # OUTPUT
%   direct: matrice contente le direzioni ottimali sulle colonne.
%   ATTENZIONE!!!: se si vuole la norma delle direzioni: vecnorm(direct)
%   teta: angolo fra una direzione ottimale e la successiva in gradi °
%
%   # INPUT
%   A : matrice
%   b : termine noto
%   x0 : guess iniziale
%   k : ultima iterata
function [direct, teta ]= grad_conj_direction(A,b,x0,max)
x = x0;
r = b - A*x;
p = r;
k = 0;
teta = [];
direct = [p];
    while k < max
        p_v = p;
        k = k + 1;
        alpha = (p'*r)/(p'*A*p);
        x = x0 + alpha*A*p;
        r = r - alpha*A*p;
        beta = ((A*p)'*r)/((A*p)'*p);
        p = r - beta*p;
        direct = [direct p];
        angolo = rad2deg(acos((p'*p_v)/(norm(p)*norm(p_v))));
        teta = [teta; angolo];
    end
