%% 2
a0 = 2500;
av = a0;
b0 = 1;
bv = b0;
n = 4;
for i = 1 : n
    an = (av + bv)/2;
    bn = sqrt(av*bv);
    av = an;
    bv = bn;

end
    Sn = 0.5*pi*a0/bn

log(4*a0)

%% 3
A = [0 3 -1;0 0 2; 7 -1 2];
Q = [0 1 0; 0 0 1; 1 0 0];
b = [1 4 5]';
R = inv(Q)*A
y = Q'*b
x = bksub(R,y)

%% 4 
%% 5
A = 17;
x = 17
for n = 1 : 6
    x = (2*x/3) + A/(3*x^2);
end
x
nthroot(A,3)
%% 6
A = [9 3 -3; 3 1 4; 2 2 1];
[L U P] = lu(A);
U
L
%% 7
%A = pentadiag ( 1, −8, 14, −8, 1 )
clear 
clc
n = 100;

A = diag(14*ones(1,n))+diag(-8*ones(1,n-1),-1)+diag(-8*ones(1,n-1),1)+diag(1*ones(1,n-2),-2)+diag(1*ones(1,n-2),2);

b = ones(n,1);
x0 = b;
[x k] = richardson(A,b,eye(n),x0,1e-3,1e5);
x(1)
k
%res = norm(r)/norm(b)
%% 
clc
clear
A = [9 -2 -2 -1; -2 7 -1 -1; -2 -1 7 -1; -1 -1 -1 5];
b = ones(4,1);
P = diag(2*ones(1,4))+diag(-ones(1,3),1)+diag(-ones(1,3),-1);
x0=zeros(4,1);
[x k] = richardson(A,b,P,x0,1e-5,2);
x = pcg(A,b,1e-3,2,P)
k
%%
S = 0;
A = 25;
N = 100;
for n = 1 : 100
    S = S + ((1/n)*((1-(1/A))^n));
end
S
%%
n = 300;
A = diag(6*ones(1,n))+diag(-2*ones(1,n-1),-1)+diag(-2*ones(1,n-1),1)+diag(ones(1,n-2),-2)+diag(ones(1,n-2),2);
Bj = eye(n)-inv(diag(diag(A)))*A;
Bgs = eye(n)-inv(tril(A))*A;
max(abs(eig(Bj)))
max(abs(eig(Bgs)))
b = A*ones(n,1);
[x k r err] = gauss_seidel(A,b,b,1e-6,1000);

k
err = norm(b-x)/norm(b)
res_norm=norm(r)/norm(b)
%%
clear 
clc
S = 0
A = 8
n = 0;
err = 1;
while err>=0.01
    n = n+1;
    S = S + (1/n)*(1-1/A)^n;
    err = log(A)-S;
end
n
S
log(A)
%%
A = [2 -1 0; -1 5 -1; 0 -1 8];
b = ones(3,1);
P = diag([4 5 6]);
alpha_opt = 2/(max(abs(eig(inv(P)*A)))+min(abs(eig(inv(P)*A))))
[x k]= richardson(A,b,P,zeros(3,1),1e-5,4,alpha_opt);
x
%%
S = 1;
n = 90;
for k = 1 :90
    S = S + ((-1)^k)/(2*k+1);
end
4*S
%%
n = 7;
b = 5*ones(n,1);
A = hilb(n);
x = A\b;
r = b-A*x;
e_rel = cond(A)*(norm(r)/norm(b))
%%
A = [10 -1 0; 0 3 5; 2 4 1];
b = [1 2 3]';
[L U P]= lu(A);
y = fwsub(L,P*b);
x = bksub(U,y);
L(2,1)
U(3,3)
y(3)
%%
n = 100;
A = diag(20*ones(1,n))+diag(-11*ones(1,n-1),-1)+diag(-11*ones(1,n-1),1)+diag(ones(1,n-2),-2)+diag(ones(1,n-2),2);
x = 2*ones(n,1);
b = A*x;
min(eig(A))
cond(A)

[x k] = richardson(A,b,eye(n),b,1e-3,1e5);
x(1)
k
res = norm(b-A*x)/norm(b)
%%
A = [7 -1 -1; -1 5 -1; -1 -1 3];
b = ones(3,1);
alpha = 2/(max(abs(eig(A)))+min(abs(eig(A))))
[x k]=richardson(A,b,eye(3),b,1e-5,4,alpha)
%%
n = 1000;
x1 = 5*ones(n,1);
A = diag(3*ones(1,n))+diag(-(3/2)*ones(1,n-1),-1)+diag(-(3/2)*ones(1,n-1),1);
b = A*x1;
[x k r]=richardson(A,b,eye(n),b,1e-2,1e3);
res_norm = norm(r)/norm(b)
err_rel = norm(x1-x)/norm(x1)
k
x(1)
[x flag res_norm k] = pcg(A,b,1e-2,1e3,[],[],b)
res_norm
err_rel = norm(x1-x)/norm(x1)
x(1)
k

%%
teta = sym('teta');
A = [5 -1 0; -1 teta 1; 0 0 3];
P = diag(diag(A));
B = eye(3)-inv(P)*A;
max(abs(eig(B)))

%%
syms k;
Sn = 4*symsum(((-1)^k)/(2*k+1),0,90);
Sn = double(Sn)

%%
A = 8;
Sn = 0;
err = 1;
n = 0;
while err>=0.01
    n = n+1;
    Sn = Sn + (1/n)*(1-1/A)^n;
    err = log(A)-Sn
end
n
%%
A = [8 -1 -1; 3 1 -1; -1 4 7];
b = [3 2 1]';
%A = sym(A);
[L U P] = lu(A);
U = sym(U)
L = sym(L)
y = fwsub(L,P*b);
y = sym(y)
%%
syms g;
A = [7 -1 -2; -1 5 -2; 0 g 1];
P = tril(A);
B = eye(3) - inv(P)*A;
abs(eig(B))

%%
syms k;
A = [2 -1; -1 -k];
pretty(eig(A))