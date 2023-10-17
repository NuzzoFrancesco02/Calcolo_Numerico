function [t,u] = Crank_Nicolson_sistemi(fun, tv, y0, Nh)
h = ( tv( end ) - tv( 1 ) ) / Nh;

u = zeros( size( y0, 1 ), Nh + 1 );
t = linspace( tv( 1 ), tv( end ), Nh + 1 );

u( :, 1 ) = y0;
for n = 1 : Nh
    phi = @(u_new) u(:,n) + h/2.*(fun(t(n),u(:,n))+fun(t(n+1),u_new));
    succ = ptofis_sys(u(:,n),phi,1e8,1e-10);
    u(:,n+1) = succ(:,end);
end