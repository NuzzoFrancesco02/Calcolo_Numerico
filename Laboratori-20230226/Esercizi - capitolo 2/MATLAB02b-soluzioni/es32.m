% clc

%% esercizio 2 - Prima parte
clear
close all

% Punto 1
N = [5:20];
K = [];
Kprec=[];
Itnp = [];
Itp = [];
for n = N;

    A = diag(4*ones(n,1)) + diag(ones(n-1, 1), 1) + diag(2*ones(n-2, 1), 2) + ...
        diag(ones(n-1, 1), -1) + diag(2*ones(n-2, 1), -2);
    KA = cond(A);
    K = [K KA];
    b = ones(n,1);
   
    toll = 1e-6;
    nmax = 5*1e3;
    P = tril(A);
    Kprec=[Kprec, cond(inv(P)*A)];

    x0 = zeros(n,1);
    
    % metodo del gradiente
    [xnp, knp] = richardson(A, b, eye(n), x0, toll, nmax);
    Itnp = [Itnp knp];
    
    % metodo del gradiente precondizionato
    [xp, kp] = richardson(A, b, P, x0, toll, nmax);
    Itp = [Itp kp];
end

% Punto 2
figure(1);
semilogy(N, Itnp, 'r', N, Itp, 'b')
grid on
title('confronto metodo del gradiente VS gradiente precondizionato')
xlabel('dimensioni sistema');
ylabel('iterazioni')
legend('gradiente', 'gradiente precondizionato', 'Location', 'Northwest')

% Punto 3
figure(2);
plot(N, K, 'r')
grid on
title('Numero di condizionamento di A')
xlabel('dimensioni sistema');
ylabel('cond(A)')

%Punto 3bis
figure(4)
plot(N, Kprec, 'r')
grid on
title('Numero di condizionamento di P^{-1}A')
xlabel('dimensioni sistema');
ylabel('cond(P^{-1}A)')


%% esercizio 2 - Seconda parte

N = [5:20];
Itcg = [];

for n = N;

    A = diag(4*ones(n,1)) + diag(ones(n-1, 1), 1) + diag(2*ones(n-2, 1), 2) + ...
        diag(ones(n-1, 1), -1) + diag(2*ones(n-2, 1), -2);
   
    b = ones(n,1);
   
    toll = 1e-6;
    nmax = 5000;
    x0 = ones(n,1);
    
    [xkcg, kcg] = conjgrad_it(A, b, x0,nmax,toll);
    Itcg = [Itcg kcg];
end

% Punto 7
figure(3);
semilogy( N, Itnp, 'r', N, Itp, 'b',N, Itcg, 'g-')
grid on
title('confronto metodo del gradiente coniugato VS gradiente e gradiente precondizionato')
xlabel('dimensioni sistema');
ylabel('iterazioni')
legend( 'gradiente','gradiente precondizionato','gradiente coniugato', 'Location', 'Northwest')