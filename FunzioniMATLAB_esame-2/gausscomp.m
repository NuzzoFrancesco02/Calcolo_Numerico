function I = gausscomp(a, b, N, f)
 
% I = gausscomp(a, b, N, f)
%
% Formula di Gauss composita a due nodi:
% Inputs:
%   a,b: estremi di integrazione,
%   N: numero di sottointervalli (m=1) formula di integrazione semplice
%   f: funzione da integrare definita come inline o anonimous
% Output:
%      I: integrale calcolato

h  = (b-a) / N;

x  = [ a : h : b ]; % vettore degli estremi dei sottointervalli

xm = [ a+h/2 : h : b ]; % vettore dei punti medi dei sottointervalli

x1 = xm - h/( 2*sqrt(3) ); % vettore contenente una parte dei 
                           % nodi di integrazione

y1 = f( x1 );

I  = h/2 * sum(y1); % prima parte dell'integrale
   
x2 = xm + h/( 2*sqrt(3) ); % vettore contenente la parte 
                           % restante dei nodi di integrazione
y2 = f( x2 );

I  = I + h/2 * sum(y2);



