%% 1

em=2^(1-6);
err = 0.5*em
%% 2
A = 8;
err = 1;
n = 0;
S = 0;
while err > 0.01
    n = n + 1;
    S = S + ((1-1/A)^n)/n;
    err = abs(log(A)-S);
end
n
%% 3
A = [8 -1 -1; 3 1 -1; -1 4 7];
b = [3 2 1]';
[L U P] = lu(A);
L = sym(L);
U = sym(U);
y = L\P*b;
L(2,1)
U(3,3)
y(3)
%% 4
A = [7 -1 -1; -1 5 -1; -1 -1 3];
b = ones(3,1);
a_opt = 2/(max(eig(A))+min(eig(A)))
x0 = b;
richardson(A,b,eye(3),x0,1e-8,4,a_opt)
%% 5
A = hilb(5);
s = 0.5;
x0 = ones(5,1);
for i = 0:2
    invpowershift(A,s,1e-8,i,x0)
end
%% 6
A = diag([3 -1 1 1])+diag([1 1 6],1);
A(4,3)=-6;
shift_finder(A,1+6i)
compass(A)
%% 7
f = @(x) log(x./6);
x = 4:0.01:8;
plot(x,f(x))
% f'(x) != 0 --> p = 2
%% 8
f = @(x) (x-6).^3;
df = Jac(f);
newton(5,1e-8,1,f,df,3)
%% 9
fatt = 0.05/(0.5^2);
k_2 = fatt*(0.05)^2
%% 10 
% 0 < mu < 1/3
%% es 1
i = 0;
lambda = @(j,n) 3 + 3.*cos(pi.*(j)./(n+1));
for n = 10:10:100
    A = diagonals([-3/2 3 -3/2],n);
    i = i + 1;
    K(i) = lambda(1,n)/lambda(n,n);
end
plot(10:10:100,K)
%% es 3
n = 1000;
A = diagonals([-3/2 3 -3/2],n);
x = 5*ones(n,1);
b = A*x;
c = rand(size(b));
c = c./norm(c);
db = c.*1e-9;
err_stim = cond(A)*norm(-db)/norm(b)

x_dx = A \ (b+db);
err_vero = norm(x_dx-x)/norm(x)
%% es 4
n = 1000;
A = diagonals([-3/2 3 -3/2],n);
x = 5*ones(n,1);
b = A*x;
[sol, iter] = richardson(A,b,eye(n),b,1e-2,1e3);
iter
sol(1)
r_norm = norm(b-A*sol)/(norm(b))
err_rel = norm(x-sol)/norm(x)
%% es 5
n = 1000;
A = diagonals([-3/2 3 -3/2],n);
x = 5*ones(n,1);
b = A*x;
[sol,~,~,iter] = pcg(A,b,1e-2,1e3,[],[],b);
sol(1)
iter
r_norm = norm(b-A*sol)/(norm(b))
err_rel = norm(x-sol)/norm(x)
%% es 6
n = 1000;
A = diagonals([-3/2 3 -3/2],n);
x = 5*ones(n,1);
b = A*x;
[p,teta] = grad_conj_direction(A,b,b,2);
teta
for i = 1 : 2
    p(:,i)'*A*p(:,3)
end
%% es 7
n = 100;
A = diagonals([-3/2 3 -3/2],n);
F = @(x) A*x + exp(-x./20) - ones(n,1);
a = zeros(n,1);
x0 = 0.1*ones(n,1);
for i = 1 : 3
    sol = broyden(F,A,x0,1e-8,i);
    sol(1)
end
%%
em=2^(1-7);
err = 0.5*em
%%
A = 216;
x = A;
for n = 1 : 10
    x = (A/(3*(x)^2)) + 2*x/3;
end
x
nthroot(A,3)
%%
A = [4 -1 -1; -1 4 -2; -1 -2 6];
b = ones(3,1);
P = diag([4 4 6]);
x0 = zeros(3,1);
a_opt = 2/(max(eig(P\A))+min(eig(P\A)))
richardson(A,b,P,x0,1e-8,5,a_opt)
%%
syms t;
A = [5 -1 0; -1 t 1; 0 0 3];
B = eye(3)-diag(diag(A))\A;
pretty(simplify(abs(eig(B))))
%%
A = hilb(7);
s = 0.2;
x0 = ones(7,1);
for i = 0:2
    invpowershift(A,s,1e-8,i,x0)
end
%%
syms t;
A = [3 t 1; -t 1 4; 0 0 9];
pretty(simplify(eig(A)))
%%
f = @(x) (cos(pi.*x)).^2;
x = 0:0.001:2;
plot(x,f(x))
%%
f = @(x) (x-1).*log(x);
df = Jac(f);
ddf = Jac(df);
newton(0.9,1e-8,1,f,df,2)
%%
fatt = 1;
k2 = fatt*(1e-2)^2;
%%
n = 1000;

A = diagonals([-1 2 -1],n);
x = ones(n,1);
b = A*x;
c = rand(size(b));
c = c./norm(c);
db = c.*1e-6;
err_stim = cond(A)*norm(-db)/norm(b)
x_dx = A\(b+db);
err_vero = norm(x_dx-x)/norm(x)
%%
n = 1000;
A = diagonals([-1 2 -1],n);
x = ones(n,1);
b = A*x;
[sol,k]= richardson(A,b,eye(n),b,1e-2,1e3);
k 
sol(1)
r_norm = norm(b-A*sol)/norm(b)
err_norm = norm(x-sol)/norm(x)
%%
n = 1000;
A = diagonals([-1 2 -1],n);
x = ones(n,1);
b = A*x;
eigpower(A,1e-6,0,ones(n,1))
eigpower(A,1e-6,1,ones(n,1))
eigpower(A,1e-6,1e3,ones(n,1))
%%
n = 1000;
A = diagonals([-1 2 -1],n);
x = ones(n,1);
b = A*x;
l_max = 2*(1+cos(pi/10001));
[l,x,it,fact] = eigpower(A,1e-6,1e3,ones(n,1),l_max);
format long;
fact
format short;
%%
n = 1000;
A = diagonals([-1 2 -1],n);
x = ones(n,1);
b = A*x;
F = @(x) exp(-x./10)+A*x-ones(n,1);
J = @(x) A-diag(exp(-x./10))./10;
x0 = ones(n,1);
for i = 1:3
    sol = gaussnewton(F,J,x0,1e-8,i);
    sol(1)
end
%%
2^(-6)
%%
a(1)=4;
b(1)=2*sqrt(2);
for i = 1:8
    a(2*i)=2*a(i)*b(i)/(a(i)+b(i));
    b(2*i)=sqrt(a(i*2)*b(i));
end
a(16)
%%
A = [9 3 -3; 3 1 4; 2 2 1];
b = [3 4 5]';
[L U P]=lu(A);
L = sym(L);
U = sym(U);
y = L\(P*b);
L(2,1)
U(3,3)
y(3)
%%
syms g;
A = [4 -1 0; -1 g 1; 0 0 1];
B = eye(3)-diag(diag(A))\A;
pretty(simplify(abs(eig(B))))
%%
A = [9 -2 -2 -1; -2 7 -1 -1; -2 -1 7 -1; -1 -1 -1 5];
b = [1 1 1 1]';
P = diagonals([-1 2 -1],4);
x0 = zeros(4,1);
x = pcg(A,b,1e-8,2,P,[],x0)
%%
d = (50-1)/(50+1);
e_k = (d^100)*1
%%
syms g;
A = [5 g; g 3];
x0 = [1 1]';
pretty(simplify(eig_approx(A,x0,0)))
%%
A = [3 -1 -1; -1 4 -1; -1 -1 5];
x0 = ones(3,1);
eigpower(A,1e-8,0,x0)
eigpower(A,1e-8,3,x0)
%%
A = [2 sqrt(10)-1 0; -sqrt(10)-1 2 0; 0 0 6];
s = 6;
shift_finder(A,s);
%%
f = @(x) sin(x+sqrt(2));
x = bisez(-2,0,1e-8,f);
x_2 = x(3)
%%

%%
A = diagonals([-1 2 -1],1000);
x = ones(1000,1);
b = A*x;
tol = 1e-2;
it = 1e3;
x0 = b;
[sol,k] = richardson(A,b,eye(1000),x0,tol,it);
k
sol(1)
res_norm=norm(b-A*sol)/norm(b)
er_rel=norm(x-sol)/norm(x)
%%
f = @(x) cos(pi.*x).*(x-0.5);
df = @(x) cos(pi.*x)-pi.*(x-0.5).*sin(pi.*x);
ddf = @(x) -pi.*sin(pi.*x)-pi.*(sin(pi.*x)+pi.*(x-0.5).*cos(pi.*x));
bis = @(x) x;
% ddf(0.5) ~= 0, mol = 2;
x0 = 0.9;
toll = 1e-6;
[x_newt, it_newt] = newton(x0,toll,1e4,f,df);
[x_newt_mod, it_newt_mod] = newton(x0,toll,1e4,f,df,2);
x_newt(end)
x_newt_mod(end)
it_newt
it_newt_mod
% newton non soddisfancente, modificato sempre soddisfacente
x_sec = secanti(f,0.9,0.7,10,1e-8);
x_sec(3)
x_sec(4)
x_sec(11)
p = stimap(x_sec);
p(end)
%%
phi = @(x,mu) x + 1.*cos(pi.*x)./(2.*pi);
dphi = @(x,mu) 1 - mu.*sin(pi.*x)./2;
bis = @(x) x;
x = -2:0.0001:2;
plot(x,bis(x),'b','LineWidth',1)
hold on;
grid on;
plot(x,phi(x,1),'r','LineWidth',2)
figure(2)
plot(x,abs(dphi(x,1)),'g','LineWidth',2)
yline(1,'LineWidth',1)
%%
F = @(x1, x2, x3) [sin(pi*x1)-x2; (x2^2)-x3; -x1-x2+(x3^2)];
J = @(x1,x2,x3) [pi.*cos(x1.*pi),-1,0;0,2.*x2,-1;-1,-1,2.*x3];
x0 = [1/5 1/5 1/5]';
newton_3var(x0,2,1e-8,F,J)
%%
A = diagonals([1 -11 20 -11 1],100);
P = diagonals([-1 2 -1],100);
K = max(eig(P\A))/min(eig(P\A))
cond(P\A)
B = P\A;
min(eig(B))
if B == B'
    disp('Sono simm')
else
    disp('Non sono simm')
end
%%
phi = @(x) x-0.5.*(exp(2.*(x-1))-1);
dphi = Jac(phi);
dphi_sol = @(x) abs(1 - exp(2.*(x-1))) -1;
bis = @(x) x;
x = -1:0.01:3;
yline(-1)
hold on;
plot(x,phi(x),'r',x,bis(x),'b')
figure(2)
yline(1)
hold on;
plot(x,abs(dphi(x)))
fsolve(dphi_sol,1.2)
%%
n = 100;
A = diagonals([0 -1 0 -2 10],n)+diag(1:n);
b = ones(100,1);
P = tril(A);
B = eye(n)-P\A;
rho_gs = max(abs(eig(B)));
toll = 1e-4;
nmax = 1000;
[sol,k] = gauss_seidel(A,b,b,toll,nmax);
k;
sol(1);
sol(2);
r_norm = norm(b-A*sol)/norm(b);
e_stim = cond(A)*r_norm
cond(A)
%%
A = hilb(3);
K = cond(A);
d = (K-1)/(K+1);
log(1/200)/log(d)
%%
A = [6 -2 -2; -2 8 -4; -2 -4 10];
b = [1 1 1]';
x0 = zeros(3,1);
pcg(A,b,1e-8,2,[],[],x0)
%%
A = [4 -3; 2 3];
x0 = [1;0];
for i = 0:1
l = eig_approx(A,x0,i)
end
%%
A = [6 0 0; 5 -1 0; 3 0 -4];
shift_finder(A,-4)
compass(eig(A))
%%
A = [6 0 0; 11 4 0; 13 7 -2];
shift_finder(A,4)
%%
f = @(x) x.*exp(x)-1;
phi = @(x) (x+1)./(exp(x)+1);
x = 0:0.01:2;
dphi=Jac(phi);
plot(x,dphi(x))
%%

A = [3 1 0 0; 0 -1 1 0; 0 0 1 6; 0 0 -6 1];
eig(A)
compass(eig(A ))
%%
phi = @(x) cos(x);
phi_d = @(x) (x*phi(phi(x))-phi(x)^2)/(phi(phi(x))+x-2*phi(x));
phi_d(1.2)
phi_d(phi_d(1.2))
%%
clear
clc
A = diagonals([-1 2 -1],1000);
x0 = ones(1000,1);
[lambda,~,~,lambda_vec] = eigpower(A,1e-6,1e3,x0);
[p,c] = stimap(lambda_vec);

%%
A = [3 2 0; 2 3 0; 0 0 8];
b = [5 5 8]';
P = diag(diag(A))
2/max(eig(P\A))
K = max(abs(eig(P\A)))/min(abs(eig(P\A)));
d = (K-1)/(K+1);
1 /(d ^ 10);
%%
A = [ 4 2 1; 2 4 1; 1 1 7];
b = [2 2 2]';
[x,~,~,alpha_0]=richardson(A,b,eye(3),b,1e-8,1);
alpha_0
x
%%
syms b;
A = [10 -2; -2 b];
x = [1;1];
pretty(simplify(eig_approx(A,x,0)))
%%
A = [4 -2; -1 1];
x0 = [1;0];
for i = 0 : 2
    eig_approx(A,x0,i)
end
%%
A = diagonals([-1 2 -1],100);
x = ones(100,1);
b = A*x;
c = rand(size(b));
c = c./norm(c);
db = c.*1e-6;
[stim, vero] = perturbated_matrix(A,x,b,db)
%%
A = [8 0 0; 1 4 0; 3 8 1];
shift_finder(A,4)
%%
f = @(x) 1-x.*exp(x);
x = -1:0.01:1;
plot(x,f(x))
grid on;
%%
phi1 = @(x,a) 0.5.*(x+a./x);
phi2 = @(x,a) x.*(3-((x).^2)./a).*0.5;
%%
n = 1000;
A = diagonals([-3/2 3 -3/2],n);
x = 5*ones(n,1);
b = A*x;
c = rand(size(b));
c = c./norm(c);
db = c*1e-9;
[stim, vero]=perturbated_matrix(A,x,b,db)
%%
A = hilb(7);
b = 5*ones(7,1);
x = A\b;
cond(A)*norm(b-A*x)/norm(b)
%%
f = @(x) exp(3*x)-1;
df = Jac(f);
ddf = Jac(df);
((1e-2)^2)*0.5*(ddf(0)/df(0))
x = newton(1e-2,1e-8,1,f,df);
abs(x(end))
%%
A = diagonals([1 -11 20 -11 1],100);

w = 1.45:0.1:1.85;
for i = 1 : length(w)
    P = tril(A)./w(i);
    B = eye(100)-P\A;
    rho(i) = max(abs(eig(B)));
end
[rho_min, pos] = min(rho);
rho_min
w(pos)
%% 
A = [6 -1 0; -1 4 -1; 0 -1 2];
x0 = ones(3,1);
eig_approx(A,x0,0)
eig_approx(A,x0,3)
%%
f = @(x) log(x./3).*sin(pi.*x);
df = Jac(f);
x = newton(2.7,1e-8,1000,f,df);
[p,c] = stimap(x);
%%
cos(cos(cos(0.9)))
%%
n = 100;
A = diagonals([1 -11 20 -11 1],n);
x = 2*ones(n,1);
b = A*x;
phi = @(y) 0.5*y'*A*y - y'*b;
phi(x)
y = ones(n,1);
norm(-b+A*y)
%%
n = 100;
A = diagonals([1 -11 20 -11 1],n);
x = 2*ones(n,1);
b = A*x;
beta = 0.5:0.001:1;
for i = 1:length(beta)
    P = diagonals([-beta(i) 2 -beta(i)],n);
    K(i) = max(eig(P\A))/min(eig(P\A));
end
[K_min, pos] = min(K);
beta(pos)
d = (K_min-1)/(K_min+1);
d^10
%%
n = 100;
A = diagonals([1 -11 20 -11 1],n);
x = 2*ones(n,1);
b = A*x;
x = b;
r = b - A*x;
p = r;
k = 0;
teta = [];
direct = [p];

    while k < 2
        k = k + 1;
        alpha = (p'*r)/(p'*A*p);
        x = x + alpha*p;
        r = r - alpha*A*p;
        beta = (p'*A*r)/(p'*A*p);
        p_v = p;
        p = r - beta*p;
        direct = [direct p];
        angolo = rad2deg(acos((p'*p_v)/(norm(p)*norm(p_v))));
        teta = [teta angolo];
    end
norm(direct(:,1))
norm(direct(:,2))
norm(direct(:,3))
teta
%%
phi = @(x) x + cos(pi.*x)./(2*pi);
dphi =  @(x) 1-sin(x.*pi)./2;
x = -2:0.01:2;
bisettrice = @(x) x;

figure(1)
plot(x,bisettrice(x),x,phi(x));
grid on
% l'intervallo del grafico che contiene 0.5 è dato dall'intersezione con la
% bisettrice :  a = -0.5  b = 1.5

figure(2)
plot(x,abs(dphi(x)))
hold on;
yline(1)
% l'intervallo in cui |dphi(x)|<1 è a = 0  b = 1 che è un intervallo più
% restrittivo
%%
f = @(x) exp(x)-2;
x = corde(f,0,2,2,2,1e-8)
%%
A = diagonals([1 -11 20 -11 1],100);
x = 2*ones(100,1);
b = A*x;
K = max(eig(A))/min(eig(A))
phi = @(y) 0.5*y'*A*y-y'*b;
phi(x)
y = ones(100,1);
norm(-b+A*y)
d = (K-1)/(K+1);
abb = 1e-3;
k_min = ceil(log(abb)/log(d))

e_0 = b-x;
e_k = (d^k_min)*sqrt(e_0'*A*e_0)
sol = richardson(A,b,eye(100),b,1e-3,1e5);
sol(1)
res_norm = norm(b-A*sol)/norm(b)
%%
A = diagonals([1 -11 20 -11 1],100);
x = 2*ones(100,1);
b = A*x;
F = @(x) A*x + exp(-2*x)-ones(100,1);
J = @(x) A + diag(-2*exp(-2*x));
x0 = 0.3*ones(100,1);

x = newton_vectorial(F,J,x0,3,1e-8);

x(1,2)
x(1,3)
x(1,4)
%%
S = 0;
tol = 0.01;
err = 1;
n = 0;
A = 8;
while err > tol
    n = n + 1;
    S = S + ((1-1/A)^n)/n;
    err = abs(S-log(A));
end
n
%%
A = [7 -1 -1; -1 5 -1; -1 -1 3];
b = ones(3,1);
a_opt = 2/(max(eig(A))+min(eig(A)))
richardson(A,b,eye(3),b,1e-8,4,a_opt)
%%
f = @(x) log(x./6);
x = 3:0.01:8;
plot(x,f(x))
%%
f = @(x) (x-6)^3;
df = Jac(f);
ddf = Jac(df);
dddf = Jac(ddf);
newton(5,1e-8,1,f,df,3)
%%
e_k = 0.5;
e_k1 = 0.05;
mu = e_k1/(e_k^2);
e_k2 = mu*(e_k1^2)
%%
n = 1000;
x = 5*ones(n,1);
A = diagonals([-3/2 3 -3/2],n);
b = A*x;
c = rand(size(b));
c = c./norm(c);
db = c*1e-9;
err_stim = cond(A)*norm(-db)/norm(b)
dx = (A \(b+db))-x;
err_es = norm(dx)/norm(x)
[sol,N] = richardson(A,b,eye(n),b,1e-2,1e3);
N
sol(1)
res_norm = norm(b-A*sol)/norm(b)
err_rel = norm(x-sol)/norm(x)
[sol,~,res_norm,N] = pcg(A,b,1e-2,1e3,[],[],b);
N
sol(1)
res_norm
err_es = norm(x-sol)/norm(x);

[dir, teta] = grad_conj_direction(A,b,b,2);
teta
for i = 1:2
    dir(:,i)'*A*dir(:,3)
end
%%
A = 216;
x = A;

for n = 1:10
    x = (A/(3*(x)^2))+2*x/3;
end
x
nthroot(A,3)
%%
A = [4 -1 -1; -1 4 -2; -1 -2 6];
P = diag([4 4 6]);
b = ones(3,1);
a_opt = 2/(max(eig(P\A))+min(eig(P\A)))
richardson(A,b,P,zeros(3,1),1e-8,5,a_opt)
%%
syms t;
A = [5 -1 0; -1 t 1; 0 0 3];
B = eye(3)-diag(diag(A))\A;
pretty(simplify(abs(eig(B))))
%%
A = hilb(7);
invpowershift(A,0.2,1e-8,2,ones(7,1))
%%
syms t;
A = [3 t 1; -t 1 4; 0 0 9];
pretty(simplify(abs(eig(A))))
%%
f = @(x) cos(pi*x)^2;
df = Jac(f);
df(0.5)
%%
f = @(x) (x-1).*log(x);
df = Jac(f);
ddf = Jac(df);
newton(0.9,1e-8,1,f,df,2)
%%
mu = (1e-2)/((1e-1)^2);
k_2 = mu*(1e-2)^2
%%
n = 1000;
A = diagonals([-1 2 -1],n);
x = ones(n,1);
b = A*x;
[sol,N] = richardson(A,b,eye(n),b,1e-2,1e3);
N
sol(1)
res_norm = norm(b-A*sol)
err_norm = norm(x-sol)/norm(x)
%%
n = 1000;
A = diagonals([-1 2 -1],n);
x0 = ones(n,1);
eigpower(A,1e-6,0,x0)
eigpower(A,1e-6,1,x0)
eigpower(A,1e-6,1e3,x0)
%%
n = 1000;
A = diagonals([-1 2 -1],n);
F = @(x) exp(-x./10)+A*x-ones(1000,1);
J = @(x) -diag(exp(-x./10))./(10) + A;
x0 = ones(1000,1);
for i = 1:3
    x = gaussnewton(F,J,x0,1e-8,i);
    x(1)
end
%%
A = 10;
err = 1;
tol = 1e-2;
S = 0;
n = 0;
while err>tol
    n = n + 1;
    S = S + ((1-1/A)^n)/n;
    err = abs(log(A)-S);
end
n
%%
n = 90;
(2*(n^3)/3) + 10 * (2*n^2)
%%
A = [9 -1 0; -1 5 -1; 0 -1 2];
b = ones(3,1);
P = diag([7 5 3]);
a_opt = 2/(max(eig(P\A))+min(eig(P\A)))
x0 = zeros(3,1);
richardson(A,b,P,x0,1e-8,4,a_opt)
%%
A = [9 4; -1 4];
qrbasic(A,1e-8,3)
%%
f = @(x) (x-3).*log(x-2).*exp(x);
df = Jac(f);
ddf = Jac(df);
x = newton(2.8,1e-8,1000,f,df,2);
stimap(x);
%%
F = @(x1,x2) [x1*x2-2;-x1+(x2^2)-2*x2+1];
x0 = [1.5;2.5];
newton_2var(x0,1,1e-8,F,Jac(F))
%%
n = 100;
A = diagonals([1 -16 30 -16 1],n);
b = ones(n,1);
P1 = diagonals([-1 4 -1],n);
P2 = diagonals([-2 4 -2],n);
K1 = max(eig(P1\A))/min(eig(P1\A));
K2 = max(eig(P2\A))/min(eig(P2\A));
d = (K2-1)/(K2+1);
abb = d^10
c = (sqrt(K2)-1)/(sqrt(K2)+1);
abb_pcg = (2*c^10)/(1+c^(2*10))

cond_spect_sp(A,P2,ones(100,1),1e-10,1)
cond_spect_sp(A,P2,ones(100,1),1e-10,2)
cond_spect_sp(A,P2,ones(100,1),1e-10,100)
%%
an = 2500;
bn = 1;
for n = 1:4
    av = an;
    bv = bn;
    an = (av+bv)/2;
    bn=sqrt(av*bv);
end
S = pi*2500/(2*bn)
%%
A = [0 3 -1; 0 0 2; 7 -1 2];
Q = [0 1 0; 0 0 1; 1 0 0];
b = [1 4 5]';
R = Q\A;
y = Q'*b;
x = R\y;
R(3,3)
y(1)
x(1)
%%
syms g;
A = [5 -1 -1; 0 4 -2; 0 g 3];
B = eye(3)-tril(A)\A;

A = [3 -g 0; g 1 0; 1 2 6];
pretty(simplify(abs(eig(A))))
%%
f = @(x) exp(x).*log((x-5)./2).*(7-x);
df = Jac(f);
df(7)
%%
phi = @(x) x*((x^2) -3*x + 3);
phi(phi(phi(1.9)))
%%
A = diagonals([1 -8 14 -8 1],100);
tol = 1e-3;
K = cond(A);
d = (K-1)/(K+1);
k = ceil(log(tol)/log(d));
b = ones(100,1);
[x,it] = richardson(A,b,eye(100),b,1e-3,1e5);
it
%%
n = 10;
A = diagonals([-1 2 -1],n);
F = @(x) A*x + sin(pi*x);
J = @(x) A +pi*diag(cos(pi*x));
x0 = (1:10)'./50;

x = newton_vectorial(F,J,x0,2,1e-10);
x(10,2)
x(10,3)
%%
A = [9 -2 -2 -1; -2 7 -1 -1; -2 -1 7 -1; -1 -1 -1 5];
b = ones(4,1);
P = diagonals([-1 2 -1],4);
x = pcg(A,b,1e-10,2,P,[],zeros(4,1))
%%
syms g;
A = [5 g; g 3];
x = [1 1]';
pretty(simplify(eig_approx(A,x,0))) 
%%
em = 1/64;
ceil(1-log(em)/log(2))
%%
S = 0;
A = 25;
for n = 1:100
    S = S + ((1-1/A)^n)/n;
end
log(A)
S
%%
n = 60;
((2*n^3)/3)+10*2*n^2
%%
b = ones(100,1);
K2 = 1e10;
r_norm = 1e-12;
e_stim = K2*r_norm/norm(b)
%%
mu = 0.005/(0.1^2);
k_2 = mu*0.005^2
%%
phi = @(x) x + (1-exp(3*x-1))/3;
dphi=Jac(phi);
ddphi=Jac(dphi);
%%
A = diagonals([1 -2 6 -2 1],300);
x = ones(300,1);
b = A*x;
[sol,k] = gauss_seidel(A,b,b ,1e-6,1000);
k
err_rel = norm(x-sol)/norm(x)
res_norm = norm(b-A*sol)/norm(b)
err_stim = cond(A)*res_norm
%%
n = 100;
A = diagonals([-1 2 -1],100);
x = ones(100,1);
b = A*x;
c = rand(size(b));
c = c./norm(c);
db = c*1e-6;
err_stim = cond(A)*norm(-db)/norm(b);
dx = (A\(b+db))-x;
err_vero = norm(dx)/norm(x);

thomas = 8*n - 7;
LU = (2*n^3)/3;

k = 1000;
e0 = zeros(100,1)-x;
K = max(eig(A))/min(eig(A))
cond(A)
d = (K-1)/(K+1);
ek_A = (d^k)*sqrt(e0'*A*e0)
%%
F = @(x) exp(x)+A*x-ones(100,1);
J = @(x) diag(exp(x)) + A;
x0 = (1:100)'./300;
x = newton_vectorial(F,J,x0,2,1e-10);
x(100,2)
x(100,3)
%%
phi = @(x) x -140.*(exp((x./7)-1)-1)./11;
x = 4:0.01:7.6672;
bi = @(x) x;
phi8 = @(x) phi(x)-7.6672;
fsolve(phi8,5)
plot(x , phi(x), x, bi(x))
figure(2)
dphi=Jac(phi);
plot(x,abs(dphi(x)));
hold on;
yline(1);
dphi1 = @(x) abs(dphi(x))-1;
fsolve(dphi1,6)
%%
f = @(x) (x-6).^3;
df = Jac(f);
ddf = Jac(df);
dddf = Jac(ddf);
m = 3;
newton(5,1e-8,1,f,df,m)
%%
f = @(x) cos(pi.*x).^2;
df = Jac(f);
x = 0:0.001:1;
plot(x,f(x))
%%
f = @(x) (x-1).*log(x);
df=Jac(f);
df(1);
ddf = Jac(df);
ddf(1);
m = 2;
newton( 0.9,1e-10,1,f,df,m)
%%
f = @(x) (x-3)*log(x-2)*exp(x);
df = Jac(f);
ddf = Jac(df);
ddf(3);
m = 2;
x = newton(2.8,1e-10,10000,f,df,2);
stimap(x);
%%
x_float(2,4,[1 0 1 1],2,0)
%%
f = @(x) exp(-4.*x)-2.*x;
x0 = 1;
it = 0;
toll = 1e-2;
err = toll+1;
x = x0;
while err > toll 
    it = it + 1;
    x = x - (f(x)^2)/(f(x+f(x))-f(x));
    err = abs(f(x));
end
it
x
%% 
t = eps_m(2,1e-8)
%%
f = @(x) exp(x)-2;
corde(f,0,2,2,2,1e-10)
%%
err = 0.5*eps_m(2,6,'eps')
%%
l = @(j,n) 3 + 3.*cos((pi.*j)./(n+1));
i = 0;
for n = 10:10:100
    i = i + 1;
    A = diagonals([-3/2 3 -3/2],n);
    K(i) = l(1,n)/l(n,n);
end
n = 10:10:100;
x = 0:1:100;
K_es = 4.*(x.^2)./(pi^2);
plot(n,K,'o',x,K_es)
%%
n = 1000;
A = diagonals([-3/2 3 -3/2],n);
x = 5*ones(n,1);
b = A*x;
[sol, it] = richardson(A,b,eye(n),b,1e-2,1e3);
it
sol(1)
res_norm = norm(b-A*sol)/norm(b)
err_norm = norm(x-sol)/norm(x)
%%
[sol,~,~,It] = pcg(A,b,1e-2,1e3,[],[],b);
It
sol(1)
res_norm = norm(b-A*sol)/norm(b)
err_norm = norm(x-sol)/norm(x)
%%
0.5*eps_m(2,7,'eps')
format shortE
%%
syms t;
A = [5 -1 0; -1 t 1; 0 0 3];
B = eye(3)-diag(diag(A))\A;
pretty(simplify(eig(B)))
%%
syms t;
A = [3 t 1; -t 1 4; 0 0 9];
pretty(simplify(abs(eig(A))))
%%
eps_m(2,0.03125,'xmin')
%%
syms g;
A = [5+4*g 10-2*g; -2*g g];
pretty(simplify(eig(A)))
%%
A = [3 -1 0; -1 2 -1; 0 -1 5];
x0 = [1 1 1]';
l_max = eigpower(A,1e-9,3,x0);
l_min = invpower(A,x0,1e-9,3);
l_max/l_min
%%
A = [11 0 0; 5 -2 0; 3 0 -9];
shift_finder(A,-9)
%%
%syms g;
g = 3;
A = [9 g 0; g 1 0; 0 0 2];

%pretty(simplify(eig(A)))
shift_finder(A,2)
%%
f = @(x) exp(-x)-1;
x = bisez(-1,4,1e-10,f);
x(1)
x(2)
%%
F = @(x1,x2) [exp(((x1^2) + (x2^2) -1)/4)-1; x1+exp(x2)-2];
J = @(x1,x2)[(x1.*exp(x1.^2./4+x2.^2./4-1./4))./2,(x2.*exp(x1.^2./4+x2.^2./4-1./4))./2;1,exp(x2)]
x0 = [2/3; 1/3];
newton_2var(x0,2,1e-10,F,J)
%%
phi = @(x) x - 140.*(exp((x./7)-1)-1)./11;
dphi = @(x) 1-(20.*exp(x./7-1))./11;

x = 4:0.01:8;
figure(1)
yline([4 8])
hold on;
grid on;
plot(x,phi(x))
% trovo il limite inferiore dove phi supera 8:
phi8 = @(x) phi(x)-8;
a = fsolve(phi8,5)

figure(2)
yline(1)
grid on;
hold on;
plot(x,abs(dphi(x)))
% trovo il limite superiore dove |dphi| supera 1
dphi1 = @(x) abs(dphi(x))-1;
b = fsolve(dphi1,7)
% noto che phi(a) supera 7.6672, quindi non è contenuta dall'intervallo:
phi7 = @(x) phi(x)-b;
a = fsolve(phi7,5)

x = a:0.01:b;
figure(3)
yline([a b])
hold on;
xline([a b])
plot(x,phi(x))
%%
A = hilb(7);
b = 5*ones(7,1);
x = A\b;
err_stim = cond(A)*norm(b-A*x)/norm(b)
%%
A = diagonals([0 -1 2 -1 0],1000);
x = 2*ones(1000,1);
b = A*x;
c = rand(size(b));
c = c./norm(c);
db = c.*1e-6;
[stim, vero] = perturbated_matrix(A,x,b,db)
