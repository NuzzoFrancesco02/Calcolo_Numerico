function [ t, u ] = teta_met(fun, tv, y0, Nh, teta)
h = ( tv( end ) - tv( 1 ) ) / Nh;

u = zeros( size( y0, 1 ), Nh + 1 );
t = linspace( tv( 1 ), tv( end ), Nh + 1 );

u( :, 1 ) = y0;
if teta == 0
    eulero_avanti_sistemi(fun, tv, y0, Nh)
end
if teta ~= 0
    for n = 1 : Nh
        phi =@(u_new) u( : , n ) + h*((1-teta)*fun(t(n),u(:,n))+teta*fun(t(n+1),u_new));
        succ = ptofis_sys(u( : , n),phi,1e10,1e-8);
        u(:, n + 1) = succ( : , end);
    end
end