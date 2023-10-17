%%
%   [x_nodes,uh] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,str1,str2,str3,h)

%
% str1 '1' : differenze finite
%      '2' : ordine due 
% str2 'D-N': Dirichelet, Neumann
%      'N-D': Neumann, Dirichelet
% ATTENZIONE: metodo differenze di ordine due per derivata prima applicato a tutti i nodi, non
% solo all'estremità

function [x_nodes,uh] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,str1,str2,str3,h)
if nargin == 12
     h = ( b - a ) / ( N + 1 );
end
Peh = abs(eta)*h/(2*mu);
if strcmp(str3,'up') 
    muh = mu*(1+Peh);
else
    muh = mu;
end
if strcmp(str2,'D-N')
   
    x_nodes = linspace( a, b, N + 2 );
    A = sparse( 1 : N, 1 : N, 2, N + 1, N + 1 ) ...
        + sparse( 2 : N, 1 : N - 1, -1, N + 1, N + 1 ) ...
        + sparse( 1 : N, 2 : N + 1, -1, N + 1, N + 1 );
    A = muh / h^2 * A;
    if strcmp(str1,'1')
    A = A + eta / ( 2 * h ) * ( sparse( 2 : N, 1 : N - 1, -1, N + 1, N + 1 ) ...
        + sparse( 1 : N, 2 : N + 1, 1, N + 1, N + 1 ) );
    elseif strcmp(str1,'2')
    A = A + eta  * ( sparse( 2 : N, 1 : N - 1, -1, N + 1, N + 1 ) ... 
        + sparse( 1 : N, 2 : N + 1, 1, N + 1, N + 1 ) ) / ( 2 * h );
    %{
    A = A + eta  * ( sparse( 2 : N, 1 : N - 1, 1, N + 1, N + 1 ) ...
        + sparse( 2 : N, 2: N, -4,N+1,N+1) +...
        + sparse( 1 : N, 2 : N + 1, 3, N + 1, N + 1 ) ) / ( 2 * h );
    %}
    end
    
    A = A + sparse( 1 : N, 1 : N, sigma( x_nodes( 2 : end - 1 ) ), N + 1, N + 1 );
    if strcmp(str1,'1')
        A( end, N : N + 1 ) = mu / h * [ -1 1 ];
    elseif strcmp(str1,'2')
        A( end, N - 1 : N + 1 ) = mu / ( 2 * h ) * [ 1 -4 3 ];
    end
    bv = ( f( x_nodes( 2 : end ) ) )';
    bv( 1 ) = bv( 1 ) + alpha * ( muh / h^2 + eta / ( 2 * h ) );
    bv( end ) = gamma;
    uh = A \ bv;
    uh = [ alpha; uh ];
elseif strcmp(str2,'N-D')

    x_nodes = linspace( a, b, N + 2 );
    A = sparse( 2 : N+1, 2 : N+1, 2, N + 1, N + 1 ) ...
        + sparse( 2 : N+1, 1 : N, -1, N + 1, N + 1 ) ...
        + sparse( 2 : N, 3 : N + 1, -1, N + 1, N + 1 );
    A = muh / h^2 * A;
    if strcmp(str1,'1')
    A = A + eta / ( 2 * h ) * ( sparse( 2 : N+1, 1 : N, -1, N + 1, N + 1 ) ...
        + sparse( 2 : N, 3 : N + 1, 1, N + 1, N + 1 ) );
    elseif strcmp(str1,'2')
    A = A + eta  * ( sparse( 2 : N+1, 1 : N, -1, N + 1, N + 1 ) ...
        + sparse( 2 : N, 3 : N + 1, 1, N + 1, N + 1 ) ) / ( 2 * h );
    %{
    A = A + eta  * ( sparse( 2 : N+1, 1 : N, 1, N + 1, N + 1 ) ...
        + sparse( 2 : N+1, 2: N+1, -4,N+1,N+1) +...
        + sparse( 2 : N, 3 : N + 1, 3, N + 1, N + 1 ) ) / ( 2 * h );
    %}
    
    end
    
    A = A + sparse( 2 : N+1, 2 : N+1, sigma( x_nodes( 2 : end -1 ) ), N + 1, N + 1 );
    if strcmp(str1,'1')
        A( 1, 1 : 2 ) = mu / h * [ -1 1 ];
    elseif strcmp(str1,'2')
        A( 1, 1 : 3 ) = mu / ( 2 * h ) * [ -3 4 -1 ];
        A( 2, 2 : 4 ) = A( 2, 2 : 4 ) + mu / ( 2 * h ) * [ -3 4 -1 ];
    end
    bv = ( f( x_nodes( 1 : end - 1) ) )';
    bv( 1 ) = gamma;
    bv( end ) = bv( end ) + alpha * ( muh / h^2 + eta / ( 2 * h ) );
    uh = A \ bv;
    uh = [ uh;alpha ];
end
if strcmp(str2,'N-N')
    x_nodes = linspace( a, b, N + 2 );
    A = sparse( 2 : N+1, 2 : N+1, 2, N + 2, N + 2 ) ...
        + sparse( 2 : N+1, 1 : N , -1, N + 2, N + 2 ) ...
        + sparse( 2 : N+1, 3 : N + 2, -1, N + 2, N + 2 );
    A = muh / h^2 * A;
    if strcmp(str1,'1')
    A = A + eta / ( 2 * h ) * ( sparse( 2 : N + 1, 1 : N, -1, N + 2, N + 2 ) ...
        + sparse( 2 : N + 1, 3 : N + 2, 1, N + 2, N + 2 ) );
    elseif strcmp(str1,'2')
    A = A + eta  * ( sparse( 2 : N + 1, 1 : N, -1, N + 2, N + 2 ) ...
        + sparse( 2 : N + 1, 3 : N + 2, 1, N + 2, N + 2 ) ) / ( 2 * h );
    %{    
        A = A + eta  * ( sparse( 2 : N + 1, 1 : N, 1, N + 2, N + 2 ) ...
        %+ sparse( 2 : N+1, 2: N+1, -4,N+2,N+2) +...
        + sparse( 2 : N + 1, 3 : N + 2, 3, N + 2, N + 2 ) ) / ( 2 * h );
    %}
    end
    if strcmp(str1,'1')
        A( 1, 1:2 ) = mu / h * [ -1 1 ];
        A( end, N + 1: N + 2) = mu / h * [ -1 1 ];
    elseif strcmp(str1,'2')
        A( 1, 1:3 ) = mu / ( 2 * h) * [ -3 4 -1 ];
        A( end, N : N + 2 ) = mu / ( 2 * h ) * [ 1 -4 3 ];
    end
    A = A + sparse( 2 : N+1, 2 : N+1, sigma( x_nodes( 2 : end - 1 ) ), N + 2, N + 2 );
    bv = ( f( x_nodes ) )';
    bv(1) = alpha;
    bv(end) = gamma;
    uh = A \ bv;
end




