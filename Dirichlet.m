%% [x_nodes,uh] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h)
% ATTENZIONE: xj = u(j+1)
%
% Peh = abs(eta)*h/(2*mu)
% mu = mu*(1+Peh)
function [x_nodes,uh] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h)
if nargin == 9
    h = (b-a)/(N+1);
end
x_nodes = linspace(a,b,N+2);
A = sparse( 1 : N, 1 : N, 2, N, N ) ...
    + sparse( 2 : N, 1 : N - 1, -1, N, N ) ...
    + sparse( 1 : N - 1, 2 : N, -1, N, N );
A = mu / h^2 * A;
A = A + (eta/(2*h))*(sparse(2:N,1:N-1,-1,N,N) + ...
    sparse(1:N-1,2:N,1,N,N));


A = A + sparse( diag( sigma( x_nodes( 2 : end - 1 ) ) ) );


bv = ( f( x_nodes( 2 : end - 1 ) ) )';
bv( 1 ) = bv( 1 ) + alpha*(mu/h^2 + eta/(2*h));
bv( end ) = bv( end ) +  beta*(mu/h^2 - eta/(2*h));
uh = A \ bv;
uh = [ alpha; uh; beta ];
