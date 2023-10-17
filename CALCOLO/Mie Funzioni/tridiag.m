function[A] = tridiag(b,a,c,n)

%Crea una matrice tridiagonale nxn con "a" principale, "b" sotto, "c" sopra

A = diag(a*ones(1,n)) + ...
	diag(b*ones(1,n-1),-1) + ...
	diag(c*ones(1,n-1),1);
end