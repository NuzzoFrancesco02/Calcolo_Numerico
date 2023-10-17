function x=bksubsimb(A,b)

% function [x] = bksub(A,b)
% Algoritmo di sostituzione all'indietro
% A: matrice quadrata triangolare superiore
% b: termine noto
% x: soluzione del sistema Ax=b

n=length(b);

if((size(A,1)~=n)||(size(A,2)~=n))
  error('ERRORE: dimensioni incompatibili')
end

if(~isequal(A,triu(A)))
  error('ERRORE: matrice non triangolare superiore')
end


if(prod(diag(A))==0)
% almeno un elemento diagonale nullo
  error('ERRORE: matrice singolare')
end

syms x

x = b(n)/A(n,n);

for i=1:n-1
   x=[(b(n-i)-A(n-i,n-i+1:n)*x)/A(n-i,n-i);x];
end
end






