function x=fwsub(A,b)

% function [x] = fwsub(A,b)
% Algoritmo di sostituzione in avanti
% A: matrice quadrata triangolare inferiore
% b: termine noto
% x: soluzione del sistema Ax=b

n=length(b);

if ((size(A,1)~=n)||(size(A,2)~=n))
  error('ERRORE: dimensioni incompatibili')
end

if ~isequal(A,tril(A))
  error('ERRORE: matrice non triangolare inferiore')
end


if (prod(diag(A))==0)
% almeno un elemento diagonale nullo
  error('ERRORE: matrice singolare')
end

x=zeros(n,1);

x(1)=b(1)/A(1,1);

for i=2:n
  x(i)=(b(i)-A(i,1:i-1)*x(1:i-1))/A(i,i);
end
