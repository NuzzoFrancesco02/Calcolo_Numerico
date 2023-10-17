clear
clc
A = [[50 1 3]; [1 6 0]; [3 0 1]];
B = [[50 1 10]; [2 20 1]; [1 4 70]];
C = [[7 8 9]; [5 4 3]; [1 2 6]];
%% 1.1.1
% controllo A
if(simm_defpos(A))
    fprintf('\nLa matrice A è simmetrica e def positiva');
end
% controllo B
if (dom_diag_stretta(B))
    fprintf('\nLa matrice B è a dominza stretta per righe o per colonne');
end
% controllo C
if(cond_nec_suff(C))
    fprintf(('\nLa matrice C rispetta la condizione necessaria e sufficiente'));
end
%% 1.1.2
% scrivo funzione laugauss.m

%% 1.1.3
% calcolo L ed U delle tre matrici
[L_A, U_A] = lugauss(A)
[L_B, U_B] = lugauss(B);
[L_C, U_C] = lugauss(C);
[L_A_es, U_A_es] = lu(A);
%% 1.1.4
% scrvo funzioni fwsub e bwsub

%% 1.1.5
% risolvo il sistema con x = [1 1 1]'
% A*x = b --> LU*x = b --> L*(U*x) = b --> L*y = b  &&  U*x = y
x = [1; 1; 1];
b = A * x;
y = fwsub(L_A, b);
y_es = fwsub(L_A_es,b);
sol = bksub(U_A, y);
sol_es = bksub(U_A_es,y);


%% 1.1.6
err_rel = norm(sol-x)/norm(x);
err_nor = norm(A*x -b)/norm(b);

%% 1.2
% A*vi = ei --> LU*vi = ei --> L*y = ei --> U*vi = y
n = size(A,1);
A1 = [];

I = eye(n);

for i = 1 : n
    y = fwsub(L_A,I(:,i));
    vi = bksub(U_A,y);
    supp = zeros(1,n);
    A1 = [A1, vi];
end

%% 3.1
clear
clc
K = 100; %N/m
L = 20; %m
N = 10;

M = -2*eye(N) + diag(ones(1,N-1),-1) + diag(ones(1,N-1),1);
term_noto = zeros(N,1);
term_noto(end) = term_noto(end)-K*L;



[L_M,L_u,sol] = thomas(M,term_noto);
sol
figure(1)
%plot(sol,'ob');
all = L*(sol(2:N)-sol(1:N-1));
%plot(sol(2:N)-sol(1:N-1),'ob');

%% 4

iter = 10;
T = zeros(5,iter);
for i = 1:iter
    N = 200*i;
    M = -2*eye(N) + diag(ones(1,N-1),-1) + diag(ones(1,N-1),1);
    b = zeros(N,1) + [zeros(N-1,1);-K*L];

    t = tic;
    x_back = M\b;
    T(1,i) = toc(t);

    t = tic;
    x_back = (-M)\(-b);
    T(2,i) = toc(t);


    t = tic;
    [LL,UU,P] = lu(M);
    x_lu = UU \ (LL \ (P*b));
    T(3,i) = toc(t);

    t = tic;
    [L_thom,U_thom,x_thom] = thomas(M,b);
    T(4,i) = toc(t);

    M_chol = -M;
    b_chol = -b;
    t = tic;
    H = chol(M_chol);
    x_chol = H\(H'\b);
    T(5,i) = toc(t);

end

hold on;
dim = 200*(1:iter);
plot(dim, T(1,:),'-ob','LineWidth',2);
plot(dim, T(2,:),'-or','LineWidth',2);
plot(dim, T(3,:),'-ok','LineWidth',2);
plot(dim, T(4,:),'-om','LineWidth',2);
plot(dim, T(5,:),'-og','LineWidth',2);
legend('backslash_non_sdp','backslash_sdp','lu','thomas','cholesky','Location','northwest')

%% 5.1
n = 20;
A = diag(4*ones(1,n)) + diag(-1*ones(1,n-1),+1) + diag(-1*ones(n-1,1),-1);
A(1,:)=1;
A(:,1)=1;
%spy(A);
[L,U,P] = lu(A);
% Per evitare il fill-in uso il pivoting totale
As = sparse(A);
[Ls,Us,Ps,Qs] = lu(As);
% L'ampiezza di banda si conserva anche in L e U (la distanza degli elementi
%non nulli dalla diagonale)

%%
K = zeros(1,10);

n = 10:10:100;
for i = 1:10
    
    A = diag(2*ones(n(i),1))+diag(-ones(n(i)-1,1),-1)+diag(-ones(n(i)-1,1),1);
    K(i)=cond(A);
end
plot(n,K)

alfa(1) = A(1,1);
e = diag(A,-1);
c = diag(A,1);
for i = 2 : n
    d(i-1) = e(i-1)/alfa(i-1);
    alfa(i) = A(i,i)-d(i-1)*c(i-1);
end
det_A = prod(alfa)
det(A)