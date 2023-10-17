%%
% u(a) = alpha
% u'(b) + b_b*u(b) = gamma
function [x_nodes,uh] = Robin(a,b,alpha,gamma,b_b,mu,eta,sigma,f,N,h,str1) 
if nargin == 10
    h = ( b - a ) / ( N + 1 );
end
    x_nodes = linspace( a, b, N + 2 );
    if strcmp(str1,'D-R')
    A = sparse( 1 : N, 1 : N, 2, N + 1, N + 1 ) ...
        + sparse( 2 : N, 1 : N - 1, -1, N + 1, N + 1 ) ...
        + sparse( 1 : N, 2 : N + 1, -1, N + 1, N + 1 );
    A = mu / h^2 * A;
    A = A + eta / ( 2 * h ) * ( sparse( 2 : N, 1 : N - 1, -1, N + 1, N + 1 ) ...
        + sparse( 1 : N, 2 : N + 1, 1, N + 1, N + 1 ) );
    A = A + sparse( 1 : N, 1 : N, sigma( x_nodes( 2 : end - 1 ) ), N + 1, N + 1 );
        A( end, N : N + 1 ) = [-mu/h mu/h+b_b];
    bv = ( f( x_nodes( 2 : end ) ) )';
    bv( 1 ) = bv( 1 ) + alpha * ( mu / h^2 + eta / ( 2 * h ) );
    bv( end ) = gamma;
    uh = A \ bv;
    uh = [ alpha; uh ];
    elseif strcmp(str1,'R-D')
    A = sparse( 1 : N, 1 : N, 2, N + 1, N + 1 ) ...
        + sparse( 2 : N, 1 : N - 1, -1, N + 1, N + 1 ) ...
        + sparse( 1 : N, 2 : N + 1, -1, N + 1, N + 1 );
    A = mu / h^2 * A; 
    A = A + eta / ( 2 * h ) * ( sparse( 2 : N, 1 : N - 1, -1, N + 1, N + 1 ) ...
        + sparse( 1 : N, 2 : N + 1, 1, N + 1, N + 1 ) );
    A = A + sparse( 1 : N, 1 : N, sigma( x_nodes( 2 : end - 1 ) ), N + 1, N + 1 );
    A( 1, 1 : 2 ) = [-mu/h mu/h+b_b];
    bv = ( f( x_nodes( 1 : end - 1 ) ) )';
    bv( 1 ) = gamma; 
    bv( end ) = bv( end ) + alpha * ( mu / h^2 - eta / ( 2 * h ) );
    uh = A \ bv;
    uh = [ uh; alpha];
    end
