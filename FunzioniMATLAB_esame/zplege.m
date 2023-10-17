function [x,w] = zplege(n,a,b)
%
% Uso:     [z,p]=zplege(n,a,b)
% Scopo:   calcola i nodi e i pesi della formula di quadratura di
%          Gauss-Legendre di ordine n
% Input:   n = ordine della formula (n+1=numero dei nodi di quadratura)
% Output:  z = vettore (colonna) contenente i nodi di quadratura
%          p = vettore (colonna) contenente i pesi della formula
%

n = n + 1;
if n <= 1
   z = 0;   p = 2;
   x = (a+b)/2 + (b-a)*z/2;
   w = (b-a)*p/2;
   return
end
jac = zeros(n);
k  = [1:n-1];
v  = k./(sqrt(4*k.^2-1));
jac = jac+diag(v,1)+diag(v,-1);
[p,z] = eig(jac);
norm2 = sqrt(diag(p'*p));    
p = (2*p(1,:)'.^2)./norm2;   
z = diag(z);

x = (a+b)/2 + (b-a)*z/2;
w = (b-a)*p/2;