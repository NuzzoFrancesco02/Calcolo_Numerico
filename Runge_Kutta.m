%%   RUNGE-KUTTA: con triu(A)==0
%
%% 
%   [t, u] = Runge_stronzo(A,b,c,fun, tv, y0, Nh)
%
%%  
%   ___________________________________
%  |                                   |
%  |         ### INPUT ###             |
%  | A,b,c : tabella di Buthcer        |
%  | fun : funzione di @(t,y)          |
%  | tv: vettore [t_0 t_max]           |
%  | y0: istante iniziale              |
%  | Nh = tf/h : numero di intervalli  |
%  |___________________________________|
%
%% 
%   _________________________________________________
%  |                                                 |
%  |              v ### OUTPUT ###                   |
%  | t : vettore che contiene gli istanti di tempo   |
%  | u : vettore che contiene le approssimazioni     |
%  |_________________________________________________|
%
%% ATTENZIONE!!! t(1) = t0 e u(1) = u0
%
%% Francesco Nuzzo
function [t, u] = Runge_Kutta(A,b,c,fun, tv, y0, Nh)
s = size(A,1);
if (~isequal(triu(A,1),zeros(s,s)))
    error('ATTENZIONE! METODO NON APPLICABILE')
end
h = ( tv( end ) - tv( 1 ) ) / Nh;
u = zeros( size( y0, 1 ), Nh + 1 );
t = linspace( tv( 1 ), tv( end ), Nh + 1 );
u( :, 1 ) = y0;

    for n = 1 : Nh
        K = [];
        for i = 1 : s
            if i == 1
                S = A(1,1);
                K_phi = @(K) fun(t(n)+c(i)*h,u(:,n)+h*S*K);
                succ = ptofis_sys(zeros(length(y0),1),K_phi,1e18,1e-12);
                K = succ(:,end);
            else
                j = find(A(i,i:end)==0,1,'first');
                if isempty(j)
                    j = s - 1;
                end
                if length(y0)==1
                    S = dot(A(i,1:i-j),K);
                else
                    S = sum(A(i,1:i-j).*K,2);
                end
                K2_phi = @(K2) fun(t(n)+c(i)*h,u(:,n)+h*(S+(A(i,end)*K2)));
                succ = ptofis_sys(zeros(length(y0),1),K2_phi,1e18,1e-12);
                K = [K succ(:,end)];
            end
        end
        u(:,n+1) = u(:,n) + h.*sum(b.*K,2);
    end
end
