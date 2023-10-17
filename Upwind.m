function [x_nodes,uh] = Upwind(a,b,alpha,gamma,mu,eta,sigma,f,N,str)

h = ( b - a ) / ( N + 1 );
Pe_h = abs(eta)*h/(2*mu);
mu = mu*(1+Pe_h);
x_nodes = linspace( a, b, N + 2 );
A = sparse( 1 : N, 1 : N, 2, N + 1, N + 1 ) ...
    + sparse( 2 : N, 1 : N - 1, -1, N + 1, N + 1 ) ...
    + sparse( 1 : N, 2 : N + 1, -1, N + 1, N + 1 );
A = mu / h^2 * A;
A = A + eta / ( 2 * h ) * ( sparse( 2 : N, 1 : N - 1, -1, N + 1, N + 1 ) ...
    + sparse( 1 : N, 2 : N + 1, 1, N + 1, N + 1 ) );
A = A + sparse( 1 : N, 1 : N, sigma( x_nodes( 2 : end - 1 ) ), N + 1, N + 1 );
if strcmp(str,'1')
    A( end, N : N + 1 ) = mu / h * [ -1 1 ];
elseif strcmp(str,'2')
    A( end, N - 1 : N + 1 ) = mu / ( 2 * h ) * [ 1 -4 3 ];
end
bv = ( f( x_nodes( 2 : end ) ) )';
bv( 1 ) = bv( 1 ) + alpha * ( mu / h^2 + eta / ( 2 * h ) );
bv( end ) = gamma;
uh = A \ bv;
uh = [ alpha; uh ];



