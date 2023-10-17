clear all
clc
clear

% esercizio 3

%% punto 1

n = 20;
I0 = 2;
R = ones(n,1);

A = -diag(R) + diag(R(1:n-1),-1);
A(1,:) = 1;

b = zeros(n,1);
b(1) = I0;

%% punto 2

[L,U,P] =lu(A);

if (P == eye(n))
    disp('non e'' stato effettuato pivoting')
else
    disp('e'' stato effettuato pivoting')
end
 
%% punto 3

y = fwsub(L, b);
i = bksub(U, y);

%% punto 4

i_ex = b(1)/n * ones(n,1);

err_rel = norm(i - i_ex)/norm(i_ex)

res_nor = norm(A*i -b)/norm(b)

Cond_A = cond(A)

%% punto 5

A_mod = A;
A_mod(2,1) = 10^3;

[L_mod,U_mod,P_mod] = lu (A_mod);

if (P_mod == eye(n))
    disp('non e'' stato effettuato pivoting')
else
    disp('e'' stato effettuato pivoting')
end

y = fwsub(L_mod, P_mod*b);
i_mod = bksub(U_mod, y);

Cond_A_mod = cond(A_mod)