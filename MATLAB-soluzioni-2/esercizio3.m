clc; clear all; close all;

x_a = 0;
x_b = 1;

p1 = @(x) x.^4 - 2 * x + 1;
p2 = @(x) 3 * x.^9 - 5 * x.^4 + 1;
g = @(x) 10 ./ ( x + 2 );
z = @(x) sqrt( x );

Iex_p1 = 1 / 5;
Iex_p2 = 3 / 10;
Iex_g = 10 * log( 3 / 2 );
Iex_z = 2 / 3;

nv = 0 : 7;
Ep1_gauss = [ ];
Ep2_gauss = [ ];
Eg_gauss = [ ];
Ez_gauss = [ ];

for n = nv
    
    % nodi di quadratura di Gauss 
    [ xi, wi ] = zplege( n, x_a, x_b );
    
    % integrali approssimati ed errori commesso
    I = sum( wi .* p1( xi ) );
    Ep1_gauss = [ Ep1_gauss, abs( I - Iex_p1 ) ];
    I = sum( wi .* p2( xi ) );
    Ep2_gauss = [ Ep2_gauss, abs( I - Iex_p2 ) ];
    I = sum( wi .* g( xi ) );
    Eg_gauss = [ Eg_gauss, abs( I - Iex_g ) ];
    I = sum( wi .* z( xi ) );
    Ez_gauss = [ Ez_gauss, abs( I - Iex_z ) ];
               
end
   
figure( 1 );
semilogy( nv, max( Ep1_gauss, 1e-15 ), 'bx-', nv, max( Ep2_gauss, 1e-15 ),'ro-', ...
          nv, max( Eg_gauss, 1e-15 ), 'ks-', nv, max( Ez_gauss, 1e-15 ),'md-', ...
          'Linewidth', 2, 'MarkerSize', 10);
xlabel('n')
ylabel('errori [log]')
grid on;
legend('Gauss, p1(x)', 'Gauss, p2(x)', ...
       'Gauss, g(x)', 'Gauss, z(x)', ...
       'Location', 'best');

%%%
Ipm = pmedcomp( x_a, x_b, 1, p1 );
Epm = abs( Ipm - Iex_p1 );
It = trapcomp( x_a, x_b, 1, p1 );
Et = abs( It - Iex_p1 );
Is = simpcomp( x_a, x_b, 1, p1 );
Es = abs( Is - Iex_p1 );

figure( 2 );
semilogy( nv, max( Ep1_gauss, 1e-15 ), 'bx-', nv, max( Epm + 0 * nv, 1e-15 ), '--', ...
          nv, max( Et + 0 * nv, 1e-15 ), '--', nv, max( Es + 0 * nv, 1e-15 ), '--', ...
          'Linewidth', 2, 'MarkerSize', 10);
xlabel('n')
ylabel('errori [log]')
grid on;
legend('Gauss, p1(x)', 'Punto Medio, p1(x)', ...
       'Trapezi, p1(x)', 'Simpson, p1(x)', ...
       'Location', 'best');

%%%
Ipm = pmedcomp( x_a, x_b, 1, p2 );
Epm = abs( Ipm - Iex_p2 );
It = trapcomp( x_a, x_b, 1, p2 );
Et = abs( It - Iex_p2 );
Is = simpcomp( x_a, x_b, 1, p2 );
Es = abs( Is - Iex_p2 );

figure( 3 );
semilogy( nv, max( Ep2_gauss, 1e-15 ), 'ro-', nv, max( Epm + 0 * nv, 1e-15 ), '--', ...
          nv, max( Et + 0 * nv, 1e-15 ), '--', nv, max( Es + 0 * nv, 1e-15 ), '--', ...
          'Linewidth', 2, 'MarkerSize', 10);
xlabel('n')
ylabel('errori [log]')
grid on;
legend('Gauss, p2(x)', 'Punto Medio, p2(x)', ...
       'Trapezi, p2(x)', 'Simpson, p2(x)', ...
       'Location', 'best');

%%%
Ipm = pmedcomp( x_a, x_b, 1, g );
Epm = abs( Ipm - Iex_g );
It = trapcomp( x_a, x_b, 1, g );
Et = abs( It - Iex_g );
Is = simpcomp( x_a, x_b, 1, g );
Es = abs( Is - Iex_g );

figure( 4 );
semilogy( nv, max( Eg_gauss, 1e-15 ), 'ks-', nv, max( Epm + 0 * nv, 1e-15 ), '--', ...
          nv, max( Et + 0 * nv, 1e-15 ), '--', nv, max( Es + 0 * nv, 1e-15 ), '--', ...
          'Linewidth', 2, 'MarkerSize', 10);
xlabel('n')
ylabel('errori [log]')
grid on;
legend('Gauss, g(x)', 'Punto Medio, g(x)', ...
       'Trapezi, g(x)', 'Simpson, g(x)', ...
       'Location', 'best');
    
%%%
Ipm = pmedcomp( x_a, x_b, 1, z );
Epm = abs( Ipm - Iex_z );
It = trapcomp( x_a, x_b, 1, z );
Et = abs( It - Iex_z );
Is = simpcomp( x_a, x_b, 1, z );
Es = abs( Is - Iex_z );

figure( 5 );
semilogy( nv, max( Ez_gauss, 1e-15 ), 'md-', nv, max( Epm + 0 * nv, 1e-15 ), '--', ...
          nv, max( Et + 0 * nv, 1e-15 ), '--', nv, max( Es + 0 * nv, 1e-15 ), '--', ...
          'Linewidth', 2, 'MarkerSize', 10);
xlabel('n')
ylabel('errori [log]')
grid on;
legend('Gauss, z(x)', 'Punto Medio, z(x)', ...
       'Trapezi, z(x)', 'Simpson, z(x)', ...
       'Location', 'best');
    