function [rho] = spectralradius(A)

% Calcola il raggio spettrale della matrice A

rho = max(abs(eig(A)));
end