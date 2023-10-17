function [ t, u ] = eulero_avanti_sistemi( fun, tv, y0, Nh )

h = ( tv( end ) - tv( 1 ) ) / Nh;

u = zeros( size( y0, 1 ), Nh + 1 );
t = linspace( tv( 1 ), tv( end ), Nh + 1 );

u( :, 1 ) = y0;

for n = 1 : Nh
    u( :, n + 1 ) = u( :, n ) + h * fun( t( n ), u( :, n ) );
        
end

return