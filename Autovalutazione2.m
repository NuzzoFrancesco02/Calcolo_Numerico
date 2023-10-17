% Autovalutazione 2 
%{
%% 1
clear
err = 1;
x = 2;

n = 0;
Sn = 1;
k = 1;

while (err >= 1e-1)
    n = n + 1;
    Sn = Sn + (x^k)/factorial(k);
    k = k + 1;
    err = exp(x)-Sn;
end

%% 2
N = 1e6;
Sn = autoval_ex2(N)

function valore = autoval_ex2(N)
    Matr = rand(N,2);
    M = 0;
    for i = 1 : N
        if (norm(Matr(i,:))<=1)
            M = M + 1;
        end
    end
    valore = 4*M/N;
end
%}
%%
M = [1 2 3; 4 5 6; 7 8 9]
[L, U,P] = lu(M);
L 
U
