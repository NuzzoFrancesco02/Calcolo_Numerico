function [T]=pentadiag(n,a,b,c,d,e)
% n dimensione della matrice
% a valore sottodiagonale -2
% b valore sottodiagonale -1
% c valore diagonale principale
% d valore sopradiaginale +1
% e valore sopradiagonale +2
T=diag(c*ones(n,1))+diag(a*ones(n-2,1),-2)+diag(b*ones(n-1,1),-1)+diag(e*ones(n-2,1),+2)+diag(d*ones(n-1,1),+1);
end