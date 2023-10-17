function [t,u] = multipasso(fun, tv, y0, Nh)
h = ( tv( end ) - tv( 1 ) ) / Nh;

u = zeros( size( y0, 1 ), Nh + 1 );
t = linspace( tv( 1 ), tv( end ), Nh + 1 );

u( :, 1 ) = y0;
u(:,2) = u(:,1)+h*fun(t(1),u(:,1));
for n = 2: Nh
    u(:,n+1)=u(:,n)+3/2*h*fun(t(n),u(:,n))-0.5*h*fun(t(n-1),u(:,n-1));
end