function [T]=tridiag(n,a,b,c)
% n dimensione della matrice
% a valore sottodiagonale
% b valore diagonale principale
% c valore sopradiagonale
T=diag(b*ones(n,1))+diag(a*ones(n-1,1),-1)+diag(c*ones(n-1,1),+1);
end