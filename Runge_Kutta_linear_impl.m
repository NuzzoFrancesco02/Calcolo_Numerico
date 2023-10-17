%%
% [t,u] = Runge_Kutta_linear_impl(RK,b,c,A,g,tv, y0, Nh)
% Risolve il caso implicito con eq. dif a coefficienti costanti
function [t,u] = Runge_Kutta_linear_impl(RK,b,c,A,g,tv, y0, Nh)
s = size(RK,1);
if (~isequal(triu(RK,1),zeros(s,s)))
    error('ATTENZIONE! METODO NON APPLICABILE')
end
h = ( tv( end ) - tv( 1 ) ) / Nh;
u = zeros( size( y0, 1 ), Nh + 1 );
t = linspace( tv( 1 ), tv( end ), Nh + 1 );
u( :, 1 ) = y0;
K = zeros(length(y0),s);
for n = 1 : Nh
    for j=1:s
        K(:,j)=(eye(length(y0)) - h*RK(j,j)*A)\(h*A*K(:,[1:j-1 j+1:end])*(RK(j,[1:j-1 j+1:end]))' + A*u(:,n) + g(t(n) + c(j)*h));        
    end
    u(:,n+1) = u(:,n) + h*sum(b.*K,2);
end
