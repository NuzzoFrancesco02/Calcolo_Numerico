function [t,u] = multipasso2(A,g, tv, y0, Nh)
h = ( tv( end ) - tv( 1 ) ) / Nh;

u = zeros( size( y0, 1 ), Nh + 1 );
t = linspace( tv( 1 ), tv( end ), Nh + 1 );

u( :, 1 ) = y0;
l = length(y0);
u(:, 2) = (eye(l)-h*A)\(u(:,1)+h*g(t(2)));
for n = 2: Nh
    u(:,n+1) = (eye(l)-2/3*h*A)\(4/3*u(:,n)-1/3*u(:,n-1)+2/3*h*g(t(n+1)));
end
