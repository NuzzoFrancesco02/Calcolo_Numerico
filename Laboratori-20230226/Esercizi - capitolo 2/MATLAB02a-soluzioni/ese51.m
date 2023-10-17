clear all
close all
clc

%% definizione della matrice A
n=100; 
A=4*diag(ones(n,1))-diag(ones(n-1,1),-1)-diag(ones(n-1,1),1);
A(1,:)=ones(1,n);
A(:,1)=ones(n,1);   

%% rappresentazione elementi non nulli di A
figure; 
spy(A); 
title('A')

%% fattorizzazione LU di A 
[L, U, P]=lu(A);

% rappresentazione elementi non nulli di L ed U
figure; 
spy(L); 
title('L')

figure; 
spy(U); 
title('U')

%% conversione delle matrici in forma sparsa
As=sparse(A); 
Ls=sparse(L); 
Us=sparse(U);

% Esaminiamo i bytes di memoria allocati
whos

%% Fattorizzazione LU con tecnica di pivoting totale (da applicare su matrici sparse)
[Ls, Us, Ps, Qs]=lu(As);

% rappresentazione elementi non nulli di L, U e P*A*Q
figure;
spy(Ls)
title ('L con pivoting totale ')

figure;
spy(Us)
title ('U con pivoting totale ')

figure;
spy(Ps*As*Qs)
title ('P*A*Q')