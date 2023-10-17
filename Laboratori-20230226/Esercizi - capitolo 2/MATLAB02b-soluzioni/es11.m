% Soluzione del Laboratorio 5 - Esercizio 1
clc
clear
close
n = 100; 
R1 = 1;
R2 = -2;


%% punto 1
A = diag(R2*ones(n,1)) + diag(R1*ones(n-1,1),-1);
A(1,:) = 1;
spy(A);
nnz(A)
Asp = sparse(A);
whos A Asp

%% punto 2

[L,U] = lugauss(A);
figure, spy(L), title('L')
figure, spy(U), title('U')
figure, spy(A), title('A')


%% punto 3

Dinv = diag(1./diag(A));
Bj = eye(n) - Dinv*A;   % matrice di iterazione di Jacobi

T = tril(A);
Bgs = eye(n) - inv(T)*A; % matrice di iterazione di Gauss-Seidel   

rho_j = max(abs(eig(Bj)))
rho_gs = max(abs(eig(Bgs)))


%% punto 5

b = ones(n,1);
b(1) = 2;
toll = 1e-6;
x0 = zeros(n,1);
nmax = 1000;
[xn,k]=jacobi(A,b,x0,toll,nmax)