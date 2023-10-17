% Autovalutazione_4
%% 1
A = [1 2; -1 4];
Pj = diag(diag(A));
Pgs = tril(A);
Bj = eye(2)-Pj\A;
Bgs = eye(2)-Pgs\A;
max(abs(eig(Bj)))
max(abs(eig(Bgs)))
%% 2 
syms g;
A = [1 g g; g 1 g; g g 1];
P = diag(diag(A));

B = eye(3)-P\A;
pretty(abs(eig(B)))
%% 3
n = 4;
b = ones(n,1);
A = diagonals([-1, 3, 2],4);
[x k] = richardson(A,b,tril(A),b,1e-6,5,1);
x
k
%% 4
syms beta;
A = [3 -1; -1 2];
b = ones(2,1);
P = [beta -1; 0 2];
B = eye(2)-P\A;
pretty(simplify(abs(eig(B))))
%% 5
A = hilb(4);
d = (cond(A)-1)/(cond(A)+1);
k = (log(1/1000))/log(d)
%% 6
clear
clc
n = 100;
A = diag(8.1*ones(n,1))+diag(-3*ones(n-1,1),-1)+diag(-3*ones(n-1,1),1)+diag(-ones(n-2,1),-2)+diag(-ones(n-2,1),2);
beta = 2:0.00001:2.5;
lambda = [];
for i = 1:length(beta)
    P = diag(beta(i)*ones(n,1))+diag(-ones(n-1,1),-1)+diag(-ones(n-1,1),1);
    alpha_opt = 2/(max(eig(P\A))+min(eig(P\A)));
    B = eye(n)-alpha_opt*(P\A);
    lambda(i) = max((eig(B)));
end
plot(beta,lambda)
