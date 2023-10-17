function rho=spectre(A)
%calcola il raggio spettrale della matrice in input
rho=max(abs(eig(A)));
end