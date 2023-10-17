function [ t, u ] = multipasso_esame( A, g, tv, y0, Nh )

h = ( tv( end ) - tv( 1 ) ) / Nh;

u = zeros( size( y0, 1 ), Nh + 1 );
t = linspace( tv( 1 ), tv( end ), Nh + 1 );

u( :, 1 ) = y0;
u(:,2) = (eye(length(y0))-h/2*A)\(u(:,1)+h/2*(A*u(:,1)+g(t(1)))+h/2*g(t(2)));
for n = 2 : Nh
    u( :, n + 1 ) = (3/2*eye(length(y0))-h*A)\(2*u(:,n)-0.5*u(:,n-1)+h*g(t(n+1)));
        
end