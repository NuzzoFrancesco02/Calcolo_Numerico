% Soluzione del Laboratorio 5 - Esercizio 2

clear, close all, clc

n = 7;

%% punto 1
A = diag(9*ones(1,n))   + diag(-3*ones(1,n-1),1) + ...
	diag(-3*ones(1,n-1),-1) + diag(ones(1,n-2),2) + ...
	diag(ones(1,n-2),-2)
b = [7 4 5 5 5 4 7]'

%% punto 2
Adiag = diag(abs(A));
Aout_diag = sum(abs(A),2) - diag(abs(A));

if (Adiag > Aout_diag) 
    disp('La matrice e'' diagonale dominante stretta per righe’) 
else
    disp('La matrice non e'' diagonale dominante') 
end

%% punto 3
if (A==A')
    if (eig(A) > 0)
        disp('La matrice e'' simmetrica e definita positiva') 
    else
        disp('La matrice e'' simmetrica ma non definita positiva')
    end
else
    disp('La matrice non e'' simmetrica') 
end

%% punto 5
toll = 1e-6;
x0 = zeros(n,1);
nmax = 1000;
[x_gs,k_gs]=gs(A,b,x0,toll,nmax)

%% punto 6
[x_jac,k_jac]=jacobi(A,b,x0,toll,nmax)

Dinv = diag(1./diag(A));
Bj = eye(n) - Dinv * A;

T = tril(A);
Bgs = eye(n) - inv(T) * A;

rho_j = max(abs(eig(Bj)))
rho_gs = max(abs(eig(Bgs)))

% stima del numero di iterazioni necessario
stima_k_j = log(toll)/log(rho_j)
stima_k_gs = log(toll)/log(rho_gs)
