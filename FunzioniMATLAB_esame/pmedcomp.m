function I = pmedcomp(a, b, N, f)

% I = pmedcomp(a, b, N, f)
%
% Formula del punto medio composita: 
% Inputs:
%      a,b: estremi di integrazione,
%      N: numero di sottointervalli (N=1 formula di integrazione semplice)
%      f: funzione da integrare definita come inline o anonimous
% Output:
%      I: integrale calcolato

h = ( b - a ) / N;

x = [ a + h / 2 : h : b - h / 2 ]; % vettore dei punti medi

I = h * sum( f(x) );


