%%
A = [4 2 -8; 2 1 8; 2 0 1];
b = [1 2 8]';
[L U P] = lu(A);
L(2,1)
U(3,3)
y = fwsub(L,P*b)
%% 
A = [8 -1 0; -1 5 -1; 0 -1 2];
b = ones(3,1);
P = diag([6 5 4]);

lambda = eig(P\A);
alpha = 2/(max(lambda)+min(lambda));
[x k] = richardson(A,b,P,zeros(3,1),1e-6,4,alpha)

%%
syms g;
A = [8 -5 g; -5 12 3; g 3 8];
solve(simplify(det(A))==0)

%%
A = [10 -2; -2 8];
x = [-1 1]';
l = invpower(A,1e-6,10,x)
%%
n = 100;
D = 1:100;
A = diag(D)+diag(-ones(n-1,1),-1)+diag(-2*ones(n-1,1),1)+diag(10*ones(n-2,1),2);
b = ones(n,1);
P = tril(A);
Bgs = eye(n)-(P\A);
max(abs(eig(Bgs)));
[x k] = gauss_seidel(A,b,b,1e-4,1000);
x(1)
x(2)
r_norm = norm(b-A*x)/(norm(b))
err_rel = cond(A)*r_norm
cond(A)
sqrt(max(eig(A'*A))/min(eig(A'*A)))
%%
syms g;
A = [2*g 2 -4; g 1 4; 2 0 1];
b = [1 2 4]';
[L U P] = lu(A);
L(2,1)
U(3,3)
y = L \ P*b;
y(2)

%%
A = [6 -2 -2; -2 8 -4; -2 -4 10];
b = ones(3,1);
P = [3 0 0; -1 4 0; -1 -2 5];
B = eye(3)-P\A;
g = P\b;
x = B*b + g
x = sym(x)
%%
A = [1 2; -1 4];
B = eye(2)-diag([1 4])\A;
abs(eig(B))
%%
A = hilb(3);
d = (cond(A)-1)/(cond(A)+1);
k = log(1/200)/log(d);
ceil(k)
%% 
A = [6 -2 -2; -2 8 -4;-2 -4 10];
b = ones(3,1);
x0 = zeros(3,1);
x = pcg(A,b,1e-6,2,[],[],x0)

%%
A = [4 -3; 2 3];
l = eigpower(A,1e-6,2,[1 0]')
%%

A = diagonals([-1 2 -1],4);
condiz = sdpcond(A,1e-8,200)
cond(A)
x0 = ones(4,1);
%x0 = [1 2 3 4]';
[V, D] = eig(A);
x0'*V(:,4);
[l_max, x] = eigpower(A,1e-8,200,ones(4,1));
%%
clear 
clc
A = diagonals([-3/2 3 -3/2],1000);
x = 5*ones(1000,1);
b = A*x;
c = rand(size(b));
c = c./norm(c);
db = c.*1e-9;
%STIMA 
err_stim = cond(A)*norm(db)/norm(b)
%VERO Aritm. esatta
x_dx = A \ (b+db);
err_ver = norm(x-x_dx)/norm(x)

A = diagonals([-3/2 3 -3/2],100);
F = @(x) A*x + exp(-x/20)-ones(size(A,1),1);
for i = 1:3
x = broyden(F,A,0.1*ones(100,1),1e-8,i);
x(1)
end
%%
A = 216;
x0 = A;
x = x0;
for n = 1 : 10
    x = A/(3*x^2)+(2/3)*x;
end
x
%%
A = [4 -1 -1; -1 4 -2;-1 -2 6];
b = ones(3,1);
P = diag([4 4 6]);
a_opt = 2/(max(eig(P\A))+min(eig(P\A)))
x = richardson(A,b,P,zeros(3,1),1e-8,5,a_opt)
%%
syms t;
A = [5 -1 0; -1 t 1; 0 0 3];
P = diag(diag(A));
B = eye(3)- P\A;
pretty(simplify(abs(eig(B))))
%%
A = hilb(7);
s = 0.2;
x0 = ones(7,1);
for i = 0:2
    l = invpowershift(A,s,1e-8,i,x0)
end
%%
syms t;
A = [3 t 1; -t 1 4; 0 0 9];
pretty(simplify(eig(A)))

%%
df = @(x) -2*pi*cos(pi*x)*sin(pi*x);
df(0.5)
%%
f = @(x) (x-1).*log(x);
df = @(x) log(x)+1-1./x;
ddf = @(x) 1./x+1./(x.^2)
x = -2:0.01:2;
%plot(x,f(x),'r',x,df(x),'b',x,ddf(x),'g');
x = newton(0.9,1e-8,1,f,df,2)

%%
n = 1000;
A = diagonals([-1 2 -1],n);
x = ones(n,1);
b = A*x;
c = rand(size(b));
c = c./norm(c);
db = 1e-6*c;
% STIMA
cond(A)*norm(db)/norm(b);
% VERO
x_cap = A\(b+db);
norm(-x+x_cap)/norm(x)
[x_cap, k, r] = richardson(A,b,eye(n),b,1e-2,1e3);
%disp('GRADIENTE')
x_cap(1);
k;
r_norm = norm(b-A*x_cap)/norm(b);
err_norm = norm(x-x_cap)/norm(x);
%disp('GRAD CONJ')
[x_cap,flag,r_norm,k]=pcg(A,b,1e-2,1e3,[],[],b);
k;
x_cap(1);
r_norm;
norm(x-x_cap)/norm(x);
l = eigpower(A,1e-6,0,ones(n,1));
l = eigpower(A,1e-6,1,ones(n,1));
[l, x ,k, fact] = eigpower(A,1e-6,1e3,ones(n,1));
fact;
%%
tol = 1e-2;
err = tol + 1;
A = 10;
N = 0;
S = 0;
while err > tol
    N = N + 1;
    S = S + ((1-1/A)^N)/N;
    err = abs(log(A)-S);
end
N
%%
A = [9 4; -1 4];
for i = 1:2
    D = qrbasic(A,1e-8,3);
    D
end
%%
i = 1;
for l = 2.1:0.001:5.9
    v(1) = 2/l;
    v(2) = l/6;
    v(3) = 6/8;
    sol(i) = max(v);
    i = i +1;
end
val = min(sol);
mini = find(sol==val);
l = 2.1:0.001:5.9;
l(mini(1))
l(mini(end))
%%
A = diagonals([1 -16 30 -16 1],100);
b = ones(100,1);
P1 = diagonals([-1 4 -1],100);
P2 = diagonals([-2 4 -2],100);
K1 = cond(P1\A)
K2 = cond(P2\A)
K2 = max(eig(P2\A))/min(eig(P2\A))
d = (K2-1)/(K2+1);
abb = d^10
c = (sqrt(K2)-1)/(sqrt(K2)+1);
(2*c^10)/(1+c^(2*10))

[x,f,r,k]=pcg(A,b,1e-3,100,[],[],b);
k
%%
A = diagonals([-1 2 -1], 100);
x = ones(100,1);
b = A*x;
%c = rand(size(b));
c = c./norm(c);
db = 1e-6*c;
stima = cond(A)*norm(db)/norm(b)
x_dx = A \ (b + db);
vero = norm(x_dx - x)/norm(x)

k = 1000;
condiz = max(eig(A))/min(eig(A))
d = (condiz-1)/(condiz+1);
x0 = zeros(100,1);
err_0 = sqrt((x-x0)'*A*(x-x0));
err_k_esimo = (d^k)*err_0

%%
n = 90;
((2*(n^3))/3) + 10*(2*n^2)

%%
a = 2500;
b = 1;
n = 0;
S = 0;
while n < 4
    a_nuovo = (a + b)/2;
    b_nuovo = sqrt(a*b);
    n = n + 1;
    a = a_nuovo;
    b = b_nuovo;

end
log(4*2500)
pi*2500/(2*b)
%% 
A = [0 3 -1; 0 0 2; 7 -1 2];
b = [1 4 5]';
Q = [0 1 0; 0 0 1; 1 0 0];
R = Q\A;
y = Q'*b;
x = R\y;
R(3,3)
y(1)
x = sym(x);
x(1)
%%
syms g;
A = [5 -1 -1; 0 4 -2; 0 g 3];
B = eye(3)-tril(A)\A;
abs(eig(B))
%%
syms g
A = [3 -g 0; g 1 0; 1 2 6];
pretty(simplify(abs(eig(A))))
%%
f = @(x) (8/(1+x^2))-x;
df = @(x) ((-16*x)/((1+x^2)^2))-1;
biseznewton(-10,10,2,2,1e-8,f,df);
%%
f = @(x) exp(x).*log((x-5)/2).*(7-x);
x = 5:0.001:9;
plot(x,f(x))
%%
phi = @(x) x .* (x.^2 -3.*x +3);
ptofis(1.9,phi,3,1e-9,1,7)
%%
A = diagonals([1 -8 14 -8 1],100);
%min(eig(A))
K = max(eig(A))/min(eig(A));
d =( K-1)/(K+1);
k = log(1e-3)/log(d);
ceil(k)
b = ones(100,1);
[x it] = richardson(A,b,eye(100),b,1e-3,1e5);
x(1)
norm(b-A*x)/norm(b)
P1 = diagonals([-1 2 -1],100);
P2 = diagonals([-2 5 -2],100);
K1 = max(eig(P1\A))/min(eig(P1\A));
K2 = max(eig(P2\A))/min(eig(P2\A));

F = @(x) exp(x) + A*x - ones(100,1);
J = @(x) A+diag(exp(x));
x0 = (1:100)'./300;
x = newton_vectorial(F,J,x0,1e-8,100);
x(end,1)
x(end,2)
%%
A = [3 -1/2 0; 1/2 1 0; 1 2 6];
x0 = ones(3,1);
for i = 0 : 2
   l = eigpower(A,1e-8,i,x0)
end
%%
A = [-4 0 0; 0 2 sqrt(10)-1; 0 -(sqrt(10)+1) 2];
%compass(eig(A))
shift_finder(A,-4);

%% 
A = [2 sqrt(17)-1 0; -(sqrt(17)+1) 2 0; 0 0 -5];
eig(A)
%gershcircles(A)
compass(eig(A))
%%
A = diagonals([-1 2 -1], 100);
F = @(x) A*x + diag(exp(x))-ones(100,1);
J = @(x) A + diag(exp(x));
x0 = [1:100]'./300;
[x, it ]= newton_vectorial(F,J,x0,100);
x(end,1)
x(end,2)
it
%% 
syms g;
A = [5 g; g 3];
x = [1,1]';
simplify(eig_approx(A,x,1))
%% 
b = ones(100,1);
err = (1e10)*(1e-12)/norm(b)
%%
A = [3 -1/2 0; 1/2 1 0; 1 2 8];
for i = 0:2
l = invpower(A,1e-8,i,ones(3,1))
end
%%

%% 
A = [7 3; -1 3];
qrbasic(A,1e-8,2)
%%
fatt = 0.005/(0.1^2)
err = (0.005^2)*fatt
%%
(30^3)/3 + 10*(2*30^2)
%%
A = [ 2 sqrt(17)-1 0; -(sqrt(17)+1) 2 0; 0 0 -5];
compass(eig(A));
shift_finder(A,-5);
%%
A = [4 2 -8; 2 1 8; 2 0 1];
b = [1 2 8]';
[L U P] = lu(A);
P
L = sym(L);
U = sym(U);
L(2,1)
U(3,3)
y = L\(P*b);
y = sym(y);
y(2)

A = [8 -1 0; -1 5 -1; 0 -1 2];
b = ones(3,1);
P = diag([6 5 4]);
l = eig(P\A);
a_opt = 2/(max(l)+min(l))
x0 = zeros(3,1);
x = richardson(A,b,P,x0,1e-8,4,a_opt)
%%
syms g;
A = [8 -5 g;-5 12 3; g 3 8];
simplify(det(A))
f = @(x) -12.*x.^2 -30.*x + 496;
x = -10:0.01:10;
plot(x,f(x))
grid on
abs(4)
%% 
A = [10 -2; -2 8];
x = [-1 1]';
eig_approx(A,x,1)
%%
A = [3 1 0 0; 0 -1 1 0; 0 0 1 6; 0 0 -6 1];
compass(eig(A))
shift_finder(A,1+6i);

%% 
iter = 0 : 90;
 4*sum (((-1).^iter)./(2.*iter + 1))
 %%
 A = 10;
 tol = 0.05;
 err = tol +1;
 n = 0;
 S = 0;
 while err > tol
    n = n + 1;
    S = (1+1/n)^n;
    err = (exp(1)-S);
 end
 n
 %%
 A = hilb(7);
 b = 5*ones(7,1);
 x = A\b;
 err = cond(A)*norm(b-A*x)/norm(b)
 %% 2
 A = [4 2 -8; 2 1 8; 2 0 1];
 b = [1 2 8]';
 [L U P] = lu(A);
 L = sym(L);
 L(2,1)
 U = sym(U);
 U(3,3)
 y = L\(P*b);
 y = sym(y);
 y(2)
 %% 3
 A = [8 -1 0; -1 5 -1; 0 -1 2];
 b = [1 1 1]';
 P = diag([6 5 4]);
 l = eig(P\A);
 a_opt = 2/(max(l)+min(l));
 x = richardson(A,b,P,zeros(3,1),1e-8,4,a_opt)
 %% 4
syms g;
A = [8 -5 g; -5 12 3; g 3 8];
pretty(simplify(det(A)));
(-15+sqrt(15^2 +12*496))/12;
(-15-sqrt(15^2 +12*496))/12;
double(solve(det(A)==0,g))
%% 5
A = [10 -2; -2 8];
x = [-1 1]';
eig_approx(A,x)
%% 6
A = [2.05 1 -1.05; 1 2 -1; -1.05 -1 2.05];
D = qrbasic(A,1e-2,100);
min(abs(D))
%% 7
% 1, zero non semplice e non è newton modificato
%% 8
phi = @(x) x - 0.5.*(exp(2.*(x-1))-1);
syms x;
x = -1:0.01:2;
plot(x,phi(x))
f = -0.5*(exp(2*(x-1))-1);

% 1.9684
%% 9
x0 = 6;
q0 = 4;
f = @(x) log((x-3)/4);
x = secanti(f,x0,q0,3,1e-8);
x
%% 
f = @(x) (x-4).^2;
x = 2:0.01:6

plot(x,f(x))
% A
%%
A = diagonals([0 -1 0 -2 10],100) + diag(1:100);
b = ones(100,1);
P = tril(A);
B = eye(100)-P\A;
% 1
rho_gs = max(abs(eig(B)))
% 2
[x,k] = gauss_seidel(A,b,b,1e-4,1000);
k
x(1)
x(2)
r_norm = norm(b-A*x)/norm(b)
% 3
err_stim = cond(A)*r_norm
% 6
[x,k] = metodoiterativo(A,b,b,1e-4,1000,0.6);
x(1)
x(2)
k
r_norm = norm(b-A*x)/norm(b)
% 7
D = diag(diag(A));
E = -tril(A,-1);
F = -triu(A,1);
w = -0.5:0.001:1.5;
i = 1;
for h = w
    B = (1-h)*eye(100) + h*(inv(D-E)*F);
    rho(i) = max(abs(eig(B)));
    i = i + 1;
end
yline(1,'-b');
hold on;
plot(w,rho)
grid on;
%%
%3 63.5
%4 28
%% 5
syms g;
A = [2*g 2 -4; g 1 4; 2 0 1];
b = [1 2 4]';
[L U P] = lu(A);
L(2,1)
U(3,3)
y = L\(P*b);
pretty(simplify(y(2)))
%% 6
A = [6 -2 -2; -2 8 -4; -2 -4 10];
b = [1 1 1]';
P = [3 0 0; -1 4 0; -1 -2 5];
x = richardson(A,b,P,b,1e-8,1,1);
x = sym(x)
%% 6
A = [6 -2 -2; -2 8 -4; -2 -4 10];
b = [1 1 1]';
P = [3 0 0; -1 4 0; -1 -2 5];
B = eye(3)-P\A;
g = P\b;
x = iterativo_generale(B,g,b,1);
x = sym(x)
%% 7
A = [1 2;-1 4];
Bj = eye(2)- diag(diag(A))\A;
Bgs = eye(2)-tril(A)\A;
rhoj = max(abs(eig(Bj)))
rhogs = max(abs(eig(Bgs)))
%% 8
A = hilb(3);
K = max(abs(eig(A)))/min(abs(eig(A)));
d = (K-1)/(K+1);
ceil(log(1/200)/log(d))
%% 9
A = [6 -2 -2; -2 8 -4; -2 -4 10];
b = [1 1 1]';
x = pcg(A,b,1e-8,2,[],[],zeros(3,1))
%% 10
A = [4 -3; 2 3];
x = [1 0]';
for i = 0:1
    eig_approx(A,x,i)   
end
    
%% 11
A = [6 0 0; 5 -1 0; 3 0 -4];
compass(eig(A))
shift_finder(A,-4);
%% 12
f = @(x) log((x-2)./3);
x = 4:0.01:6;
plot(x,f(x))
x = corde(f,4,6,4,2,1e-8)
%%
f = @(x) (x-2).*exp(x-1);
phi = @(x) x - ((x-2).*exp(x-2))./(2*exp(1)-1);
dphi = @(x) 1 - (exp(x-2).*(x-1))./(2*exp(1)-1);
x = 1.5:0.001:2.5;
plot(x,phi(x),'b',x,(dphi(x)),'r')
hold on;
yline(1)
dphi(2)
%x = ptofis(1.5,phi,1000,1e-4,1,3);
%x(2)
%x(3)
phi_delta = @(x) (x.*phi(phi(x))-(phi(x).^2))./(phi(phi(x))+x-2.*phi(x));
%x = ptofis(1.5,phi_delta,1000,1e-4,1,3);
%x(2)
%x(3)
%p = stimap(x)
%%
syms g;
A = [1 g g; g 1 g; g g 1];
B = eye(3) - diag(diag(A))\A;
abs(eig(B))
%%
A = diagonals([1 4 1],10);
l = eig(A);
a_opt = 2/(max(l)+min(l));
B = eye(10)-a_opt*A;
max(abs(eig(B)))
%%
A = [5 -1 -1; -1 7 -1; -1 -1 9];
b = [1 1 1]';
x = pcg(A,b,1e-8,2,[],[],zeros(3,1))
%%
A = [6 0 0; 11 4 0; 13 7 -2];
compass(eig(A))
shift_finder(A,4);
%%
A = [3 7 -1; 2 0 -1; 0 -1 4];
D = qrbasic(A,1e-2,100)
compass(eig(A))
%%
A = diag([-6i -4i -2i 0 2i 4i 6i]);

compass(eig(A))
l = invpowershift(A,7i)
shift_finder(A,4i)
%%
A = diagonals([-1 2 -1],100);
s = 2 + 2*cos(pi*47/101);
compass(eig(A))
shift_finder(A,s)
%%
f = @(x) exp(x).*log((x-5)./2).*(7-x);
x = 5:0.01:8;
plot(x,f(x))
grid on;
%%
syms e x;
A = [2 -2 0; (e-2) 2 0; 0 -1 3];
[L U P] = lu(A);
limit(L(3,2),e,0,'right')
%%
syms g;
A = [2*g 2 -8; g 1 8; 2 0 1];
b = [1 2 8]';
[L U P] = lu(A);
L(2,1)
U(3,3)
y = L \ (P*b);
y(2)
%%
A = [6 -2 -2; -2 8 -4; -2 -4 10];
b = ones(3,1);
P = [6 0 0; -1 8 0; -1 -2 10];
B = eye(3)-P\A;
g = P\b;
iterativo_generale(B,g,b,2)
%%
A = [3 2 0; 2 3 0; 0 0 8];
b = [5 5 8]';
a_max1 = 2/(max(eig(A)))
P = diag(diag(A))
a_max2 = 2/max(eig(P\A))
max(abs(eig(eye(3)-diag(diag(A))\A)))
K = max(eig(P\A))/min(eig(P\A));
d = (K-1)/(K+1);
d^10
%%
A = [4 2 1; 2 4 1; 1 1 7];
b = [2 2 2]';
[x,~,~,alpha]=richardson(A,b,eye(3),b,1e-8,1);
alpha
%%
syms b;
A = [10 -2; -2 b];
x = ones(2,1);
simplify(eig_approx(A,x))
%%
A = [4 -2;-1 1];
x0 = [1 0]';
for i = 0:2
eigpower(A,1e-8,i,x0)
end
%% 
A = [8 0 0; 1 4 0; 3 8 1];
shift_finder(A,4)
invpowershift(A,4)
%%
f = @(x)1 -x.*exp(x);
x = -2:0.01:2;
plot(x,f(x))
grid on
%%
A = [-4 0 0; 0 2 sqrt(10)-1; 0 -sqrt(10)-1 2];
eig(A)
%%
phi1 = @(x) 0.5.*(x+a./x);
phi2 = @(x) (x./2).*(3-((x.^2)./a));
dphi1 = @(x) 0.5-0.5*a./(x.^2);
%dphi2 = @(x) 3/2-(3/(2*a)).*()
%%
A = [10 -2; -2 1];
eig(A)
%%
A = [4 -2; -1 1];
x0 = [1, 0]';
for i = 0: 2
invpower(A,x0,1e-8,i)
end
%%
syms g,
assume(g,'real');
A = [10 0 0; 1 g 0; 7 0 2];
shift_finder(A,g)
%%
syms g;
A = [3 -g 0; g 1 0; 1 2 6];
det(A)
pretty(simplify(eig(A)))
%%
A = diagonals([1 -11 20 -11 1], 100);
Bj = eye(100)-diag(diag(A))\A;
Bgs = eye(100)-tril(A)\A;
max(abs(eig(Bj)));
max(abs(eig(Bgs)));
b = 5*ones(100,1);
tol = 1e-2;
nmax = 1e4;
[x,k,r] = gauss_seidel(A,b,b,tol,nmax);
k;
x(1);
r_norm = norm(r)/norm(b);
cond(A)
err_stim = cond(A)*r_norm
P1 = eye(100);
P2 = diagonals([-1 2 -1],100);
K1 = max(abs(eig(P1\A)))/min(abs(eig(P1\A)))
K2 = max(abs(eig(P2\A)))/min(abs(eig(P2\A)))
d1 = (K1-1)/(K1+1);
d2 = (K2-1)/(K2+1);

%abb1 = d1^10
abb2 = d2^10

c = (sqrt(K2)-1)/(sqrt(K2)+1)
abb_gj = (2*c^(10))/(1+c^(2*10))
%%
f = @(x) exp(x)-1;
a = -1;
b = 1;
for i = 1 : 2
corde(f,a,b,1,i,1e-8)
end
%%
n = 100;
A = diagonals([-2 5 3], n);
b = ones(n,1);
Pj = diag(diag(A));
Bj = eye(n) - Pj\A;
rho_j = max(abs(eig(Bj)));
x0 = b;
nmax = 1000;
tol = 1e-3;
[x, k] = jacobi(A,b,x0,tol,nmax);
k
x(1)
x(2)
r_norm = norm(b-A*x)/norm(b)
%%
f = @(x) exp(-4.*x)-2.*x;
x = 1;
k = 0;
err = 1;
while err > 1e-2 
    x = x - (f(x)^2)/(f(x+f(x))-f(x));
    k = k + 1;
    err = abs(f(x));
end
x
k
plot(-1:0.01:1,f(-1:0.01:1))
%% 
phi = @(x) x - 140.*(exp((x./7)-1)-1)./11;
x = 4:0.01:8;
dphi = @(x) 1 - 20.*(exp((x./7)-1))./11;
plot(x,phi(x),'r',x,abs(dphi(x)),'g')
yline([1 4 8])
grid on;
%%
A = diagonals([1 -11 20 -11 1],100);
b = 0.5:0.0001:1;
j = 1;
% for i = linspace(0.5,1,5000)
%     P = diagonals([-i 2 -i],100);
%     K(j) = max(abs(eig(P\A)))/min(abs(eig(P\A)));
%     j = j + 1;
% end
[~,i] = min(K);
b(i)
x = 2*ones(100,1);
b = A*x;

[p,teta] = grad_conj_direction(A,b,b,2);
vecnorm(p)
teta
%%
S = 0;
n = 0;
err = 1;
A = 8;
while err > 0.01
    n = n + 1;
    S = S + ((1-1/A)^n)/n;
    err = abs(log(A)-S);
end
n
%%
A = [3 1 0 0; 0 -1 1 0; 0 0 1 6; 0 0 -6 1];
eig(A)
shift_finder(A,1+6i)
%%
f = @(x) log(x./6)
x = 3:0.01:10;
plot(x,f(x))
%% 
n = 1000;
A = diagonals([-3/2 3 -3/2],n);
x = 5 * ones(n,1);
b = A*x;
[~,teta] = grad_conj_direction(A,b,b,2);
teta
%%
f = @(x) (x-1).*log(x);
df = @(x) log(x)+1-1./x;
ddf = @(x) (1./x) + (1./(x^2));
newton(0.9,1e-8,1,f,df,2)

F = @(x1,x2) [x1.*x2-2; -x1+(x2.^2)-2.*x2+1];
dF = @(x1,x2) []
x0 = [1.5, 2.5]';
newton(x0,1e-8,1,F)
%%
F = @(x1,x2) [exp(((x1.^2)+(x2.^2)-1)./4)-1;x1+exp(x2)-2];
J = Jac(F,'[x1,x2]');
x0 = [2/3,1/3]';
newton_2var(x0,2,1e-8,F,J)
f = @(x) 3*x^2 + 2*x + 5;
%%
f = @(x) (1-exp(x-3)).*(((x.^2)./9)-1);
df = Jac(f);
ddf = Jac(df);
x = 0:0.01:6;
subplot(2,1,1)
plot(x,f(x),'r',x,df(x),'b')
grid on
subplot(2,1,2)
plot(x,f(x),'r',x,ddf(x),'g')
grid on
%%
%xmin=B^(L-1);
B = 2;
xmin = 0.03125;
L = (log(xmin)/log(B))+1
%%
x = Floating_Point(2,4,-5,5,[1 0  1 1],0,2)
%%
B = 2;
em = 1e-8;
t = ceil(1 - log(em)/log(B))
%%
[~,~,~,~,~,~,~,err] = Floating_Point(2,6,-7,7)
%%
A = diagonals([1 -8 14 -8 1],100);
K = max(eig(A))/min(eig(A));
d = (K-1)/(K+1);
k = ceil(log(1e-3)./log(0.9998))
%%
syms g;
A = [3, 3; g 1];
x0 = [0, 1]';
l = eig_approx(A,x0,1);
pretty(simplify(l))
%%
phi = @(x) x -140.*(exp((x./7)-1)-1)./11;
dphi = @(x) 1 - 20.*(exp((x./7)-1))./11;
x = 4:0.01:8;
plot(x,phi(x),'r',x,abs(dphi(x)),'g')
grid on;
hold on;
yline([1 4 8]);
fsol
%%
n = 0:90;
S = 4*sum(((-1).^n)./((2.*n)+1))
%%
A = diagonals([-1 2 -1],100);
shift_finder(A,2.2173)
%%
A = diagonals([1 -11 20 -11 1],100);
Bj = eye(100)-diag(diag(A))\A;
Bgs = eye(100)-tril(A)\A;
rho_j = max(abs(eig(Bj)));
rho_gs = max(abs(eig(Bgs)));
b = 5 * ones(100,1);
x = gauss_seidel(A,b,b,1e-2,1e4);
err_rel = cond(A)*norm(b-A*x)/norm(b);
cond(A);
T = tril(A);
w = 1.45:0.1:1.85;
k = 1;
for i = w
    Pw = (1/i)*T;
    Bw = eye(100)-Pw\A;
    %g = Pw\b;
    rho(k) = max(abs(eig(Bw)));
    k = k + 1;
end
min(rho);
rho;
P1 = eye(100);
P2 = diagonals([-1 2 -1],100);
K1 = max(eig(P1\A))/min(eig(P1\A));
K2 = max(eig(P2\A))/min(eig(P2\A));
if isequal(P2',P2) && min(eig(P2))>0
    disp('simm e def pos')
end
d = (K2-1)/(K2+1);
k = 10;
abb = d^k;
c = (sqrt(K2)-1)/(sqrt(K2)+1);
abb = (2*c^(k))/(1+c^(2*k))
%%
A = hilb(7);
s = 0.2;
x0 = ones(7,1);
for i = 0:2
    invpowershift(A,s,1e-8,i,x0)
end
%%
A = [10 -1 0; 0 3 5; 2 4 1];
b = [1 2 3]';
[L, U, P] = lu(A);
L = sym(L);
U = sym(U);
y = L\P*b;
x = U\y;
L(2,1)
U(3,3)
y(3)
%%
A = [6 -1 0; -1 4 -1; 0 -1 2];
x0 = ones(3,1);
eigpower(A,1e-8,0,x0)
eigpower(A,1e-8,3,x0)
%%
A = [-4 0 0; 0 2 sqrt(10)-1; 0 -sqrt(10)-1 2];
shift_finder(A,-4)
%%
f = @(x) log(x./3).*sin(pi.*x);
df = Jac(f);
x = 1:0.01:5;
plot(x,f(x),'r',x,df(x),'g')
grid on;
hold on;