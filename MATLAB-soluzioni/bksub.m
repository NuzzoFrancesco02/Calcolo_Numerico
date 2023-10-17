function x=bksub(A,b)

% function [x] = bksub(A,b)
% Algoritmo di sostituzione all'indietro
% A: matrice quadrata triangolare superiore
% b: termine noto
% x: soluzione del sistema Ax=b

x=[];
n=length(b);

if((size(A,1)~=n)||(size(A,2)~=n))
  error('ERRORE: dimensioni incompatibili')
  return
end

if(A~=triu(A))
  error('ERRORE: matrice non triangolare superiore')
  return
end


if(prod(diag(A))==0)
% almeno un elemento diagonale nullo
  error('ERRORE: matrice singolare')
  return
end

x=zeros(n,1);

x(n) = b(n)/A(n,n);

for i=n-1:-1:1
   x(i)=(b(i)-A(i,i+1:n)*x(i+1:n))/A(i,i);
end

