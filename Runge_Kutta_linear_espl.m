%%
% [t,u] = Runge_Kutta_linear_espl(RK,b,c,A,g,tv, y0, Nh)
% Risolve il caso eslpicito con eq. dif a coefficienti costanti
function [t,u] = Runge_Kutta_linear_espl(RK,b,c,A,g,tv, y0, Nh)
s = size(RK,1);
if (~isequal(triu(RK,1),zeros(s,s)))
    error('ATTENZIONE! METODO NON APPLICABILE')
end
h = ( tv( end ) - tv( 1 ) ) / Nh;
u = zeros( size( y0, 1 ), Nh + 1 );
t = linspace( tv( 1 ), tv( end ), Nh + 1 );
u( :, 1 ) = y0;
for n = 1 : Nh
    K = zeros(length(y0),s);
    for i = 1 : s
        S = 0;
        for j = 1 : s
            if j<i
                S = S + RK(i,j)*K(:,j);
            end
        end
        K(:,i) = g(t(n)+c(i)*h)+A*(u(:,n)+h*S);
    end
    u (:,n+1) = u(:,n)+h*sum(b.*K,2);
end