function [ t, u ] = Giulio( A, g, tv, y0, Nh )
% Metodo di eulero all'indietro per sistemi con funzione in forma:
% f = A * y(t) + g(t)

h = ( tv( end ) - tv( 1 ) ) / Nh;

u = zeros( size( y0, 1 ), Nh + 1 );
t = linspace( tv( 1 ), tv( end ), Nh + 1 );

u( :, 1 ) = y0;
[L,U] = lu( eye(length(y0)) - h*A );
for n = 1 : Nh
    x=fwsub(L,(u(:,n) + h*g(t(n + 1))));
    u( :, n + 1 )=bksub(U,x);
end

return