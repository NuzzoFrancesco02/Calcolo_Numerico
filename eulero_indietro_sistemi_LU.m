function [ t, u ] = eulero_indietro_sistemi_LU(A, g, tv, y0, Nh)
h = ( tv( end ) - tv( 1 ) ) / Nh;

u = zeros( size( y0, 1 ), Nh + 1 );
t = linspace( tv( 1 ), tv( end ), Nh + 1 );

u( :, 1 ) = y0;
M = (eye(length(y0))-h*A);

[L,U,P] = lu(M);
for n = 1:Nh
    a = fwsub(L,P*(u(:,n)+h*g(t(n+1))));
    u(:,n+1) = bksub(U,a);
end