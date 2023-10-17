%Autovalutazione_3
%%
n = 100;
A = diag(2*ones(1,n))+diag(-ones(1,n-1),1);
oke = ones(n,1);
A(n,1)=1;
A = sparse(A);
[L U P Q] = lu(A);
%Ax = b
%P*A*Q = P*L*U*Q 
%P*L*U*Q*Q^(-1)*x=P*b
%L_*U_*Q^(-1)*x = P*b
%
y = fwsub(L,P*oke);
x_ = bksub(U,y);
L(2,1)
y(n)
x_(n)

% Ax=b ---> P*A*Q = P*L*U*Q ---> P*A*Q*Q^(-1)*x=P*b 
% P*L*U*Q*Q^(-1)*x = P*b 
% L_*U_*Q^(-1)*x = P*b  x_ = Q^(-1)*x
% L_*y_ = P*b
% U*x_ = y_
% x = Q*x_ 
%%
syms g;
A = [2*g sqrt(3)*g/2; sqrt(3)*g/2 g];
eig(A)

%%
A = [6 -2 -2; -2 8 -4; -2 -4 10];
oke = [1 1 1]';
P = diag([6 8 10]) + diag([-1 -2],-1) + diag(-1,-2);
    x0=oke;
    toll = 1e-6;
    nmax=2;
    K = 0;
    r = oke-A*x0;
    err = norm(r)/norm(oke);
    x = x0;
    B = eye(3)-(inv(P))*A;
    g = inv(P)*oke;
    while err > toll && K < nmax
        K = K+1;
        x = B*x + g;
        r = oke-A*x;
        err = norm(r)/norm(oke);
    end
x

%%
clear
clc
A = [3 1 1; 4 6 5; 8 8 9];
oke = [1 2 3]';
x0 = [1 1 1]';
toll = 1e-9;
nmax = 1000;
[x, K] = jacobi(A,oke,x0,toll,nmax);
x
xp=A\oke

%%
A = [6 -2 -2; -2 8 -4; -2 -4 10];
oke = ones(3,1);
P = [6 0 0; -1 8 0; -1 -2 10];
B = eye(3)-inv(P)*A;
g = inv(P)*oke;
err = 1;
K = 0;
x = oke;
while err > 1e-3 && K < 2
    K = K + 1;
    x = B*x + g;
    err = norm(oke-A*x)/norm(oke);
end
x

%%
A = [1 2;-1 4];
disp('jacobi');
B = eye(2)-inv(diag(diag(A)))*A;
max(abs(eig(B)))
disp('G-S');
B = eye(2)-inv(tril(A))*A;
max(abs(eig(B)))
%%
syms g;
A = [1 g g; g 1 g; g g 1];
B = eye(3)-inv(diag(diag(A)))*A;
abs(eig(B))
%%
n = 4;
A = diag(3*ones(1,n))+diag(-ones(1,n-1),-1)+diag(2*ones(1,n-1),1);
oke = ones(n,1);
x0 = oke
[x K] = gauss_seidel(A, oke, x0, 1e-5, 5)
%%
syms oke;
A = [3 -1; -1 2];
P = [oke -1; 0 2];
B = eye(2)-inv(P)*A;
abs(eig(B))
%%
A = hilb(4);
P = eye(4);
K = cond(inv(P)*A);
d = (K-1)/(K+1);
k = -3/(log10(d))
%%
clc
clear
n = 10;
A = diag(8.1*ones(1,n))+diag(-3*ones(1,n-1),-1)+diag(-3*ones(1,n-1),+1)...
    +diag(-ones(1,n-2),-2)+diag(-ones(1,n-2),2);
syms b;
P = diag(b*ones(1,n))+diag(-ones(1,n-1),-1)+diag(-ones(1,n-1),1);
B = eye(n)-inv(P)*A;
max(abs(eig(B)))


