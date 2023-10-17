function [x_nodes,uh] = trasporto_diffusione(a,b,alpha,beta,mu,eta,f,N)
h = ( b - a ) / ( N + 1 );
x_nodes = linspace( a, b, N + 2 );
A = sparse( 1 : N, 1 : N, 2, N, N ) ...
    + sparse( 2 : N, 1 : N - 1, -1, N, N ) ...
    + sparse( 1 : N - 1, 2 : N, -1, N, N );
A = mu / h^2 * A;
A = A + eta / ( 2 * h ) * ( sparse( 2 : N, 1 : N - 1, -1, N, N ) ...
    + sparse( 1 : N - 1, 2 : N, 1, N, N ) );
bv = ( f( x_nodes( 2 : end - 1 ) ) )';
bv( 1 ) = bv( 1 ) + alpha * ( mu / h^2 + eta / ( 2 * h ) );
bv( end ) = bv( end ) + beta * ( mu / h^2 - eta / ( 2 * h ) );
uh = A \ bv;
uh = [ alpha; uh; beta ];