%% 13
f = @(x) cos(3*x+sqrt(5));
a = 0; b = 1; n = 4;
poly_stima(f,n,a,b,'equi')
%% 14
f = @(x) exp(x);
a = -1; b = 2; H = 1; N = (b-a)/H;
interp_tratti(f,a,b,2,N,'equi',1.5)
%% 15
x = 0:10;
f = @(x) exp(x./10)+0.1.*sin(pi.*x+sqrt(2));
polyval(polyfit(x,f(x),2),11)
%% 16
% C
%% 17
f = @(x) exp(x);
a = 0; b = 4; 
gausslegendre_comp(a,b,f,1,2,'equi')
%% 18
syms g x;
f = g*x^3-4;
df = jacobian(f,x);
ddf = jacobian(df,x);
dddf = jacobian(ddf,x);
h = 1/3;
e = -1/12*h^2*2*dddf
%% 19
f = @(t,y) 3*exp(y)-19*t; y0 = 0;
h = 0.1; tf = 10; tv = [0 tf];
[t,u] = eulero_avanti(f,10,y0,h);
u(2)
%% 20
mu = 1; eta = 0; sigma = @(x) 0.*x + 5; f = @(x) 1 + 0.*x;
a = 0; b = 1; alpha = 1; beta = 0;
h = 0.1; N = 9;
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N);
u(6)
%% 21
e1 = 8*1e-3; h1 = 0.1; h2 = h1/2;
e2 = e1*(h2/h1)^2
%% 22
a = 0; b = 1;
mu = 1; f = @(x,t) 0; u_s = @(t) 0; u_d = @(t) 0; g_0 = @(x) 3*sin(pi*x); T = 10; h = 1/2;
delta_t = 1/8; theta = 1;
[u, x, t] = EqCalore_DiffFin_Theta(mu, f, a, b, u_s, u_d, g_0, T, h, delta_t,theta);
u(2,2)
%% ESERCIZIO
clear
clc

f = @(t) 10;
A = [-2 -6; 1 0];
g = @(t) [f(t) 0]';
fun = @(t,y) A*y + g(t);
tf = 5; tv = [0 tf];  y0 = [1 4]'; h = 0.1;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(2,2),u(2,end)
x_es = @(t) exp(-t).*(7/3.*cos(sqrt(5).*t)+10/(3*sqrt(5)).*sin(sqrt(5).*t))+5/3;
h = [10 5 2.5 1.25].*1e-3;
Eh = [];
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh abs(u(2,end)-x_es(tf))];
end
Eh,stimap_2(Eh,h);
loglog(h,Eh,h,h,'-',h,h.^2,':','LineWidth',2);
legend('Eh','H','H^2');
h = 0.1;
[t,u] = Heun_sistemi(fun,tv,y0,tf/h);
u(2,2),u(2,end)
l = max(eig(A));
options = optimoptions('fsolve','Display','off');
fsolve(@(h) abs(1+h*l+(h*l)^2/2)-1,1,options)
h = 0:0.001:1;
figure()
plot(h,abs(1+h.*l+(h.*l).^2./2),'LineWidth',2), yline(1)
%% 13
f = @(x) x+sin(pi.*x+1); n = 10; a = 0; b = 5;
x_nod = a:(b-a)/n:b;
polyval(polyfit(x_nod,f(x_nod),1),6)
%% 14
f = @(x) exp(x); a = 0; b = 2; e = 1e-2;
M = ceil((b-a)*sqrt(exp(b)/(8*e)))
x_dis = a:0.001:b; x_nod = a:(b-a)/M:b;
y_dis = interp1(x_nod,f(x_nod),x_dis);
max(abs(y_dis-f(x_dis)))
%% 15
x = 1:5; y = [1 1 0 1 2];
spline(x,y,4.5)
%% 16
syms beta g;
f = @(x) beta*x.^3+g*x.^2+1; a = -2; b = 2; M = 2;
pmedcomp(a,b,M,f)
%% 17 
a = -1; b = 1; h = 1/2;
n = (b-a)/h-2;
r = n + 1
%% 18
f = @(x) 2.^x; h = 1/4;
df = @(x) (-3*f(x)+4*f(x+h)-f(x+2*h))/(2*h);
df(0)
%% 19
f = @(t,y) 1 + t*y; y0 = 1; tf = 10; tv = [0 tf];
h = 0.1;
[t,u] = Crank_Nicolson(f,tf,y0,h);
u(2)
%% 20
mu = 1; eta = 0; sigma = @(x) 1+10.*x; f = @(x) 0.*x + 5;
a = 0; b = 1; alpha = 0; beta = 0;
h = 0.1; N = 9;
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N);
u(6)
%% 21
h1 = 1e-1; h2 = 1e-2; e1 = 1e-3;
e2 = e1*(h2/h1)^2
%% 22
mu = 0.5; f = @(x,t) 0;
a = 0; b = 1; u_s = @(t) 0; u_d = @(t) 0; g_0 = @(x) 5*sin(pi.*x);
delta_t = 0.1; h = 0.5; theta = 0; T = 10*delta_t;
[u, x, t] = EqCalore_DiffFin_Theta(mu, f, a, b, u_s, u_d, g_0, T, h, delta_t,theta);
u(2,end)
%% ESERCIZIO
m = 1; c = 2; k = @(t) 3*(2-exp(-t)); 
z = @(t) 10.*exp(-t).*cos(t).*(30.*exp(-t).*(2-exp(-t)).*cos(t)-2);
fun = @(t,y) [-c/m.*y(1) - k(t)./m.*y(2).^2 + z(t)./m; y(1)];
y0 = [-10 10]'; tf = 2; tv = [0 tf];
h = 0.05;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(:,end)
A = @(t,y) [-c/m -k(t)/m*y(2); 1 0];
g = @(t) [z(t)/m;0];
fun = @(t,y) A(t,y)*y+g(t);
[t,u] = provaboh(A,g,tv,y0,tf/h);
u(:,end)
y_es = @(t) [-10.*exp(-t).*(cos(t)+sin(t)); 10*exp(-t).*cos(t)];
Eh = [];
h = [10 5 2.5 1.25].*1e-4;
for i = h
    [t,u] = provaboh(A,g,tv,y0,tf/i);   
    Eh = [Eh norm(u(:,end)-y_es(tf))];
end
Eh, stimap_2(Eh,h);
x_es = @(t) 10*exp(-t).*cos(t);
J = @(t) [-c/m -2*k(t)/m*x_es(t); 1 0];
l = [];
for i = 0:0.001:tf
    l = [l max(eig(J(i)))];
end
l = max(l);
fsolve(@(h) abs(1+h*l+(h*l)^2/2)-1,1)
%%
clear
clc
h = 0.5;
fun = @(u) -1/h^2*(2*(0-u)-1*(u));
fsolve(@(h) fun(h)-5,1)
%%
f = @(x) cos(3.*x + sqrt(5));
n = 4; a = 0; b = 1; poly_stima(f,n,a,b,'equi')
%%
f = @(x) exp(x); a = -2; b = 2; H = 1;
interp_tratti(f,a,b,2,4,'equi',1.5)
%%
x = 0:10;
f = @(x) exp(x./10)+0.1.*sin(pi.*x+sqrt(2));
polyval(polyfit(x,f(x),2),11)
%%
x = linspace(0,1,50);rng(1);
y = 2*x.^2+0.2*sin(100*pi*x)+0.2*randn(1,50);
polyval(polyfit(x,y,2),0.5)
%%
a = 0; b = 10; f = @(x) (x+2).^2; 
x_nod = a:(b-a)/5:b;
interp1(x_nod,f(x_nod),1)
%%
f = @(x) sin(10*x)-7*x;
a = 0; b = 3; 
df = Jac(f); ddf = Jac(df);
x_dis = a:0.001:b;
df_max = max(ddf(x_dis)); tol = 1e-4;
M = ceil((b-a)*sqrt(df_max/(tol*8)))
%%
a = -1; b = 1;
n = 6; 
[eq, phi, A, w] = poly_equation(n,a,b,'CGL');
w = matlabFunction(w);
pmedcomp(a,b,1,w)
%%
f = @(x) exp(x);
integr_stima(-2,2,1e-1,f,'t')
%%
A = [0 0; 3/4 0]; b = [1/3 2/3]; c = [0 3/4];
fun = @(t,y) -2*(t+1)*y^2; tf = 10; y0 = 3;
tv = [0 tf]; h = 0.1;
[t,u] = Runge_Kutta(A,b,c,fun,tv,y0,tf/h);
u(2)
%%
f = @(t,y) -5*(t+1)*y; y0 = 8; tf = 100; tv = [0 tf];
h = 0.2; [t,u] = Crank_Nicolson(f,tf,y0,h);
u(2)
%%
A = [-3 -4; 1 0]; y0 = [7 1]';
l = max(eig(A));
fsolve(@(h) abs(1+h*l)-1,1)
%%
mu = 3; eta = -50;
h = 2*mu/abs(eta)
%%
K1 = 1e4; h1 = 1; h2 = h1/10;
K2 = K1*(h1/h2)^2
%%
h = 0.1;
N = 9;

alpha = 5;
beta = 0;
a = 0;
b = 1;
mu = 1; eta = 100; f = @(x) 0.*x;
% Peh = abs(eta)*h/(2*mu);
%mu = mu*(1+Peh);
xnodes = linspace(a, b, N+2);
Peh = (abs(eta) * h) / (2 * mu);
muUW = mu*(1+Peh);
AD = diagonals([-1 2 -1],N);
AD = (muUW / (h^2)) * AD;
AT = diagonals([-1 0 1],N);
AT = (eta/(2*h)) * AT;
A = AD + AT;
bv = f((xnodes(2:end-1)'));

bv(1) = bv(1) + (muUW/h^2+eta/(2*h)) * alpha;

bv(end) = bv(end) + (muUW/h^2-eta/(2*h)) * beta;
uh = A \ bv;
uh = [alpha;uh;beta];
uh(10)
%%
clear
clc
h = 0.1;
N = 9;

alfa = 5;
beta = 0;
a = 0;
b = 1;

eta = 100;
mu = 1;

f = @(x) 0 * x + 0;

xnodes = linspace(a, b, N+2);

Peh = (abs(eta) * h) / (2 * mu)

muUW = mu * (1 + Peh)

% AA = sparse(1:N, 1:N, 2, N, N) + sparse(2:N, 1:N-1, -1, N, N) + sparse(1:N-1,...
%     2:N, -1, N, N);

AA = diagonals([-1 2 -1],N);

AB = (muUW / (h^2)) * AA;

% BA = sparse(2:N, 1: N-1, -1, N, N) + sparse(1:N-1, 2:N, 1, N, N);

BA = diagonals([-1 0 1],N);

BB = (eta/(2*h)) * BA;

A = AB + BB;

bv = f((xnodes(2:end-1)'));
bv(1) = bv(1) + (((muUW/(h^2))+(eta/(2*h))) * alfa);
bv(end) = bv(end) + (((muUW/(h^2))-(eta/(2*h))) * beta);

uh = A \ bv;

uh = [alfa; uh; beta]
uh(end-1)
%%
h = 0.1;
N = 9;

alfa = 5;
beta = 0;
a = 0;
b = 1;
eta = 100;
mu = 1;
f = @(x) 0 .* x ;

xnodes = linspace(a, b, N+2);

Peh = (abs(eta) * h) / (2 * mu)

muUW = mu * (1 + Peh)

% AA = sparse(1:N, 1:N, 2, N, N) + sparse(2:N, 1:N-1, -1, N, N) + sparse(1:N-1,...
%     2:N, -1, N, N);

AA = diagonals([-1 2 -1],N);

AB = (muUW / (h^2)) * AA;

% BA = sparse(2:N, 1: N-1, -1, N, N) + sparse(1:N-1, 2:N, 1, N, N);

BA = diagonals([-1 0 1],N);

BB = (eta/(2*h)) .* BA;

A = AB + BB;

bv = f((xnodes(2:end-1)'));
bv(1) = bv(1) + (muUW/h^2+eta/(2*h)) * alfa;

bv(end) = bv(end) + (((mu/(h^2))-(eta/(2*h))) * beta);
uh = A \ bv;
uh = [alfa; uh; beta]
uh(end-1)
%%
f = @(x) 5+x.^7;
gausscomp(0,3,3,f)
gausslegendre_comp(0,3,f,3,1,'equi')