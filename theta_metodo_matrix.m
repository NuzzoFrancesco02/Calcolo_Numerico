%% FUNZIONE THETA-METODO

% f: sistema di funzioni
% tv : [t0 tf]
% h: paso h
% y0: dato iniziale
% theta: numero da 0 a 1

% Attenzione la funzione va bene solo se la matrice A presenta dei termini
% scalari ovvero la funzione è a coefficenti costanti

function [t,u] = theta_metodo_matrix(A,g,tv,h,y0,theta)

Nh = (tv(end)-tv(1))/h;
t = linspace(tv(1),tv(end),Nh+1);
u = zeros( size( y0, 1 ), Nh + 1 );
u( :, 1 ) = y0;
A_theta = sparse(eye(size(A,1))-h*theta*A);

for it=2:Nh+1
    
   
    bv = (eye(size(A,1))+h*(1-theta)*A)*u(:,it-1)+h*((1-theta)*g(t(it-1))+theta*g(it));
    u(:,it) = A_theta \ bv;

end