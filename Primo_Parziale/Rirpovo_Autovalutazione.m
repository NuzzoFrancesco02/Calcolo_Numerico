%% RIPROVO TEST AUTOVALUTAZIONE

%%
clc
clear
Sn = 1;
toll = 1e-1;
x = 2;
err = 1;
k = 0;
while err > toll
    k = k + 1;
    Sn = Sn + (x^k)/factorial(k);
    err = exp(2)-Sn;
end
k
%%
N = 1e6;
num = rand(2,N);
M = sum(vecnorm(num)<=1);
S = 4*M/N

%%
clear 
clc
syms g
A = [2*g 2 -8; g 1 8; 2 0 1];
b = [1 2 8]';
[L U P] = lu(A);
y = L\(P*b);
%y = fwsub(L,P*b);
% x = bksub(U,y);
L(2,1)
U(3,3)
y(2)
%%
syms g;
A = [2*g sqrt(3)*g/2; sqrt(3)*g/2 g];
abs(eig(A))
%%
clear 
clc
A = [6 -2 -2; -2 8 -4; -2 -4 10];
b = [1 1 1]';
P = [6 0 0; -1 8 0; -1 -2 10];
k = 0;

B = eye(3)-inv(P)*A;
g = inv(P)*b;
x = b;
while k < 2 
    k = k + 1;
    x = B*x + g;
end
x

%%
n = 90;
k = 0 : n;
S = 4 * sum(((-1).^k)./(2.*k +1))

%%
n = 80;
10*(n^2)

%%
A = [10 -1 0; 0 3 5; 2 4 1];
b = [1 2 3]';
[L U P] = lu(A);
L(2,1)
U(3,3)
y = fwsub(L,P*b);
y(3)
x = bksub(U,y);
A*x
%%
n = 100;
A = diag(20*ones(n,1))+diag(-11*ones(n-1,1),-1)+diag(-11*ones(n-1,1),1)+diag(ones(n-2,1),-2)+diag(ones(n-2,1),2);
x = 2*ones(n,1);
b = A * x;
min(eig(A));
K = max(eig(A))/min(eig(A));

phi = @(x) 0.5*x'*A*x - x'*b;
phi(x)
norm(A*ones(n,1)-b)

%%
n = 100;
A = 20*eye(n) + diag(ones(n-2,1),-2) ...
- 11*diag(ones(n-1,1),-1) - 11*diag(ones(n-1,1) ,1) ...
+ diag(ones(n-2,1),2);
b = 5*ones(n,1);
T = tril(A);

%%
n = 100;
A = diag(20*ones(n,1))+diag(-11*ones(n-1,1),-1)+diag(-11*ones(n-1,1),1)+diag(ones(n-2,1),-2)+diag(ones(n-2,1),2);
D = 20*diag(ones(1,n));
Bj = eye(n)-inv(D)*A;
rho_j = max(abs(eig(Bj)));