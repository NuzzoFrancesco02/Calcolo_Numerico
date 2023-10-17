function I = trapcomp(a, b, N, f)

% I = trapcomp(a, b, N, f)
%
% Formula del punto medio composita: 
% Inputs:
% Formula del trapezio composita: 
%    a,b: estremi di integrazione,
%    N: numero di sottointervalli (m=1 formula di integrazione semplice)
%    f: funzione da integrare definita come inline o anonimous
% Output:
%    I: integrale calcolato

h = ( b - a ) / N;

x = [ a : h : b ]; % vettore dei nodi di integrazione

y = f( x );
   
I = h * ( 0.5 * y( 1 ) + sum( y( 2 : N ) ) + 0.5 * y( N + 1 ) );


