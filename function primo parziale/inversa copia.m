function Ainv=inversa(A)
%restituisce la matrice inversa di A
n=size(A,1);
E=diag(ones(n,1));
X=lusolve(A,E);
Ainv=X';
end