%%
x = 0:0.5:1;
y = [2 1 1.5];
P = polyfit(x,y,2);
dP = polyder(P);
polyval(P,0.25)
polyval(dP,0.25)
%%
x = 0:0.25:1;
y = [3 0.5 1.5 -0.5 1];
x_dis = 0:0.0001:1;
min(interp1(x,y,x_dis))
%%
f = @(x) sin(pi.*x);
interp_tratti(f,0,2,1,4,'equi',1.6)
%%
x = [0 0.25 0.5 0.75 1];
y = [3 0.5 1.5 -0.5 1];
poly_scarto_IMQ(x,y,2)
%%
f = @(x) sqrt(2+abs(x));
N = 3;
trapcomp(-1,2,N,f)
%%
M1 = 20;
e1 = 1e-1; e2 = 1e-3;
M2 = ceil(M1*nthroot(e1/e2,2))
%%
f = @(t,y) -(1+sin(t))*y^2/81;
df = @(t,y) -2*(1+sin(t))*y/81;
tf = 5; tv = [0 tf]; y0 = 9; h = 0.2;
[t,u] = eulero_indietro_newton(f,df,tf,y0,h);
u(end)
%%
mu = 1; sigma = @(x) 0.*x + 2; f = @(x) 3 + x; eta = 0;
a = 0; b = 1; alpha = 1; gamma = 0;
h = 0.1; N = 9;
[x,u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,'2','D-N');
plot(x,u)
%%
mu = 1; eta = -50; sigma = @(x) 0.*x; f = sigma;
a = 0; b = 1; alpha = 0; beta = 3;
Peh = abs(eta)*h/(2*mu);
mu = mu*(1+Peh);
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,9,0.1);
plot(x,u)
u(2)
%%
mu = 3;
%%
K = 1;
A = @(y) [-1 -(147/16+K*y(2)); 1 0];
g = @(t) [-9*exp(-t/2)*(sin(3*t)^2-1)-9/2*exp(-t/4)*sin(3*t);0];
tf = 5; tv = [0 tf]; y0 = [-3/4 3]';
fun = @(t,y) A(y)*y +g(t);
h = 1e-2;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(:,2), u(:,end)
[u_min,i] = min(u(2,:));
u_min, t(i)

y = @(t) 3*exp(-t/4).*[-1/4*cos(3*t)-3*sin(3*t) cos(3*t)]';
h = [10 5 2.5 1.25].*1e-4;
Eh = [];
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh norm(u(:,end)-y(tf))];
end
Eh
stimap_2(Eh,h);

h = 1e-2;
A = [0 0; 2/3 0]; b = [1/4 3/4]; c = [0 2/3];
[t,u] = Runge_Kutta(A,b,c,fun,tv,y0,tf/h);
u(:,2), u(:,end)

A =  [-1 -(147/16); 1 0];
l = max(eig(A));
R_abs = @(h) abs(1 + h*l + (h*l)^2/2);
fsolve(@(h) R_abs(h)-1,1)
fun = @(t,y) A*y;
[t,u] = Crank_Nicolson_sistemi(fun,tv,y0,tf/h);
u(:,2),u(:,end)
%% 1
f = @(x) sin(x+sqrt(2));
[y_pol, P, err] = poly_lagr(f,3,0,pi,2,'equi');
y_pol
dP = polyder(P);
polyval(dP,2)
%%
f = @(x) exp(3*x); poly_stima(f,3,-1,1,'CGL')
%%
f = @(x) ((x/pi).^2).*(x>=0 & x<pi)+((2-x/pi).^2).*(x>=pi & x<2*pi);
i = 0:5; x = 2*pi/5.*i;
spline(x,f(x),x(4))
%%
x = 0:4; y = [2 2 0 1 0];
poly_minim_IMQ(1,x,y)
%%
f = @(x) exp(x); a = -1; b = 3; N = (b-a)/1;
trapcomp(a,b,N,f)
%%
f = @(t,y) -2*y;
y0 = 7; teta = 3/4; h = 0.1; tf = 10;
[t,u] = teta_met(f,[0 tf],y0,tf/h,teta);
u(2)
%%
e1 = 4*1e-3; h1 = 1e-2; e2 = 1e-3;
h2 = sqrt(e2/e1)*h1
%%
mu = 1; eta = 0; sigma = @(x) 1 + 19*x; f = @(x) 0.*x + 5;
N = 9; h = 0.1;
[x,u] = miste(0,1,0,0,mu,eta,sigma,f,N,'1','N-N');
u(6)
%%
mu = 1; eta = 0; sigma = @(x) 0.*x; f = @(x) 0.*x + 6;
a = 0; b = 1; b_b = 3; gamma = 0; N = 9; h = 0.1;
[x,u] = Robin(a,b,alpha,gamma,b_b,mu,eta,sigma,f,N,h);
u(end)
%%
f = @(y) [-2*y(1)-10*y(2)^2 y(1)]';
z = @(t) exp(-t/2)*(2*cos(t)-7/2*sin(t))+40*exp(-t)*sin(t)^2;
g = @(t) [z(t) 0]';
fun = @(t,y) f(y)+g(t);

tf = 10; tv = [0 tf]; y0 = [2 0]'; h = 1e-2;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(2,2),u(2,end)
plot(t,abs(u(2,:)),'LineWidth',2); yline(0.2,'LineWidth',2)
i = find(abs(u(2,:))>=0.1,1,'last'); t(i+1)

x = @(t) 2*exp(-t./2).*sin(t);
h = [10 5 2.5 1.25].*1e-4;
Eh = [];
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh max(abs(u(2,:)-x(t)))];
end
Eh, stimap_2(Eh,h);
A = @(t) [-2 1; -20*x(t) 0];
l = [];
t = 0:0.001:tf;
for i = 1:length(t)
    l = [l max(eig(A(t(i))))];
end
l = max(l);
R_abs = @(h) abs(1+h.*l);
fsolve(@(h) R_abs(h)-1,1)
[t,u] = multipasso(fun,tv,y0,tf/1e-2);
u(2,2),u(2,3),u(2,end)
%%
f = @(x) 1 + x.^5 + abs(x);
trapcomp(-1,1,1,f)
%%
V = 1;
mu = 1; eta = V; sigma = @(x) 0.*x; f = @(x) sin(x) + cos(x);
a = 0; b = 1; alpha = 0; beta = cos(1);
h = 0.1; N = (b-a)/h-1;
x_nodes = linspace(a,b,N+2);
A = zeros(N+1);
A = A + diagonals([-1 2 -1],N+1).*mu./(h^2);
A = A + eta*diagonals([-1 0 1],N+1)./(2*h);
A(end,:)=zeros(1,N+1);
A(end,end-1:end) = eta*[-1 1]./(h);
A = sparse(A);
bv = f(x_nodes(2:end)');
bv(1) = bv(1) + alpha*(mu/h^2 + eta/(2*h));
bv(end) = beta;
u = A\bv;
u = [alpha;u];
u(end)
max(abs(u-sin(x_nodes)'))
%%
P = polyfit(x_nodes,u,4);
f = @(x) P(1).*x.^4+P(2).*x.^3+P(3).*x.^2 + P(4).*x + P(5);
trapcomp(0,1,1e5,f)
h = [0.1 0.05 0.025 0.00125];
u_ex = @(x) sin(x);
errv = [];
gamma = beta;
f = @(x) sin(x) + cos(x);
Eh = [];
for i = h
N = (b-a)/i-1;
x_nodes = linspace(a,b,N+2);
A = zeros(N+1);
A = A + diagonals([-1 2 -1],N+1).*mu./(i^2);
A = A + eta*diagonals([-1 0 1],N+1)./(2*i);
A(end,:)=zeros(1,N+1);
A(end,end-1:end) = eta*[-1 1]./(i);
A = sparse(A);
bv = f(x_nodes(2:end)');
bv(1) = bv(1) + alpha*(mu/i^2 + eta/(2*i));
bv(end) = beta;
u = A\bv;
u = [alpha;u];
Eh = [Eh max(abs(u-sin(x_nodes)'))];
end
Eh,log(Eh(end)/Eh(end-1))/log(h(end)/h(end-1))
%%
clear
clc
V = 100;
mu = 1; eta = V; sigma = @(x) 0.*x; f = @(x) 0.*x;
a = 0; b = 1; alpha = 0; beta = cos(1);
h = 0.1; N = (b-a)/h-1;
x_nodes = linspace(a,b,N+2);
A = zeros(N+1);
A = A + diagonals([-1 2 -1],N+1).*mu./(h^2);
A = A + eta*diagonals([-1 0 1],N+1)./(2*h);
A(end,:)=zeros(1,N+1);
A(end,end-1:end) = eta*[-1 1]./(h);
A = sparse(A);
bv = f(x_nodes(2:end)');
bv(1) = bv(1) + alpha*(mu/h^2 + eta/(2*h));
bv(end) = beta;
u = A\bv;
u = [alpha;u];
plot(x_nodes,u)
Peh = abs(eta)*h/(2*mu);
mu = mu*(1+Peh);
A = zeros(N+1);
A = A + diagonals([-1 2 -1],N+1).*mu./(h^2);
A = A + eta*diagonals([-1 0 1],N+1)./(2*h);
A(end,:)=zeros(1,N+1);
A(end,end-1:end) = eta*[-1 1]./(h);
A = sparse(A);
bv = f(x_nodes(2:end)');
bv(1) = bv(1) + alpha*(mu/h^2 + eta/(2*h));
bv(end) = beta;
u = A\bv;
u = [alpha;u];
plot(x_nodes,u)
u(end)
%%
f = @(x) x.^2 + abs(x);
n = 5; h  = ( 2/ 4); x_nod = -1:h:1;
x_dis = -1:0.001:1;
y_dis = polyval(polyfit(x_nod,f(x_nod),2),x_dis);
trapz(x_dis,y_dis)
%%
T = 0.5; P1 = [0 0]; P2 = [1 0]; P3 = [0 1];
P = 1/3*(P1+P2+P3);
f = @(x) exp(2*x(1)+x(2));
I = abs(T)*f(P)
%%
syms t h y;
f = -t*y;
EAI_sym(f,1,5)
%%
V = 1;
a = 0; b = 1; alpha = 1; gamma = exp(1);
mu = 1; eta = V; sigma = @(x) 0.*x + 1; f = @(x) exp(x);
h = 0.1; N = (b-a)/h-1;
[x,u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,'1','D-N');
u(2),u(end)
trapz(x,u)
h = [0.1 0.05 0.025 0.0125];
Eh = [];
u_ex = @(x) exp(x);
for i = h
    N = (b-a)/i-1;
    [x,u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,'1','D-N');
    Eh = [Eh max(abs(u'-u_ex(x)))];
end
Eh
stimap_2(Eh,h);
h = 0.1; N = (b-a)/h-1;
V = 1000; eta = V; sigma = @(x) 0.*x; f = @(x) 0.*x;
[x,u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,'1','D-N',h);
plot(x,u)

%%
clear
clc

a = 0; b = 1; alpha = 1; gamma = exp(1); eta = 1000;
sigma = @(x) 0.*x; f = sigma;
mu = 1; h = 0.1; N = (b-a)/h-1;

[x,u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,'1','D-N','up');
u(end)
%%
clear
clc
f = @(x) exp(x).*(1+sin(pi.*x)); a = -1; b = 1;
n = 5; h = (b-a)/n;
x_nod = a:h:b; 
x_dis = a:0.001:b;
y_dis = polyval(polyfit(x_nod,f(x_nod),n),x_dis);
eh = max(abs(y_dis-f(x_dis)))
[~,~,err] = poly_lagr(f,5,a,b,a:0.001:b,'equi')
%%
T = [6 5.5 6 4.5 4 3.5 4 4.5 5];
day = 1:length(T);
polyval(polyfit(day,T,2),10)
%%
syms g; M = 10; a = 0; b = 1;
%e = H^2/8*max(f'')

f = @(x) x^g;
df = Jac(f,'[x]');
ddf = @(x) (g-1).*g.*x.^(g-2);
H = 1/M;
e = H^2/8*ddf(1)
%%
f = @(x) sqrt(1+x);
trapcomp(-1,3,4,f)
%%
M1 = 10; e1 = 1e-1; M2 = 100;
e2 = e1*(M1/M2)^4
%%
f = @(x,y) exp(x+3*y);

eps = [-1/sqrt(3) 1/sqrt(3)];
a = 0; b = 1; c = 0; d = 1;
x = (a+b)/2+(b-a)/2.*eps;
y = (c+d)/2+(d-c)/2.*eps;
S = 0;
for i = 1:2
    for j = 1:2
        S = S + f(x(i),y(j));
    end
end
I = (b-a)*(d-c)*S/4
%%
syms h;
f = @(t,y) -(1+t)*y; y0 = 2;
u_star = y0 + h*f(0,y0);
u = y0 + 0.5*h*(f(0,y0)+f(h,u_star));
expand(u)
%%
mu = 1; sigma = @(x) 0.*x; eta = 0; f = @(x) (2+x).^2;
a = 0; b = 1; alpha = 1; beta = 0; N = 9; 
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N);
u(6)
%%
mu = 1; eta = 40; sigma = @(x) 0.*x; f = sigma;
a = 0; b = 1; alpha = 7; beta = 0;
h = 0.1; N = 9;
Peh = abs(eta)*h/(2*mu);
mu = mu*(1+Peh);
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
u(10)
%%
A = diagonals([3 -2 -1],9);
l = max(eig(A));
fsolve(@(h) abs(1 + h*l)-1,1)
g = @(t) exp(sin(pi*t))*ones(9,1); y0 = 4*ones(9,1);
fun = @(t,y) A*y + g(t);
tf = 10; tv = [0 tf]; h = 0.1;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(5,2)
u(5,end)
[u_min,i] = min(u(5,:));
u_min, t(i)
%%


y5 = 1.142174435142178;
h = [10 5 2.5 1.25].*1e-3;
Eh = [];
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh abs(u(5,end)-y5)];
end
Eh, stimap_2(Eh,h);
h = 0.1;
[t,u] = Crank_Nicolson_sistemi(fun,tv,y0,tf/h);
u(5,2)

A = [1/3 0; 1/3 1/3]; b = [1/2 1/2]; c = [1/3 2/3];
[t,u] = Runge_Kutta(A,b,c,fun,tv,y0,tf/h);
u(5,2), u(5,end)
%%
x = linspace(0,1,50);
rng(1);
y = 2*x.^2+0.2*sin(100*pi*x)+0.2*randn(1,50);
polyval(polyfit(x,y,2),0.5)
%%
a = 0;  b = 10; f = @(x) (x+2).^2;
n = 5; h = (b-a)/n;
x_nodes = a:h:b;
interp1(x_nodes,f(x_nodes),1)
%%
f = @(x) sin(10*x)-7*x;
a = 0; b = 3;
df = Jac(f); ddf = Jac(df);
n = ceil(sqrt(3^2/8*max(ddf(a:0.001:b))/1e-4))
%%
clear
clc
a = -1; b = 1;
n = 6;
syms x;
n = 6;
w = 1;
k = 0:n;
t = -cos(pi*k/n);
x_nod = ((b-a)/2)*t + (a+b)/2;
for i = 1 : n + 1
    w = w * (x-x_nod(i));
end
w = matlabFunction(expand(w));
pmedcomp(-1,1,1,w)
%%
f = @(x) exp(x);
a = -2; b = 2; tol = 1e-1;
integr_stima(a,b,tol,f,'t')
%%
clear
clc
A = [0 0; 3/4 0]; b = [1/3 2/3]; c = [0 3/4];
fun = @(t,y) -2*(t+1)*y^2;
tf = 10;
tv = [0 tf]; y0 = 3; h = 0.1;
[t,u] = Runge_Kutta(A,b,c,fun,tv,y0,tf/h);
u(2)
%%
fun = @(t,y) -5*(t+1)*y;
y0 =  8; tf = 100;
tv = [0 tf]; h = 0.2;
[t,u] = Crank_Nicolson(fun,100,y0,h);
u(2)
%%
A = [-3 -4; 1 0];
fun = @(t,y) A*y; y0 = [7 1]';
l = max(eig(A));
fsolve(@(h) abs(1+h*l)-1,1)
%%
mu = 3; eta = -50; sigma = @(x) 0.*x; f = @(x) 0.*x;
a = 0; b = 1; alpha =1; beta = 0;
h = 2*mu/abs(eta)
%%
eta = 0;
a = 0; b = 1; mu = 2; sigma = @(x) 0.*x +1;
f = @(x) 1-x+(4+2*pi^2)*cos(pi/2*x);
alpha = 0; gamma = -2;
h = 0.1;
[x,u]= miste(a,b,alpha,gamma,mu,0,sigma,f,(b-a)/h-1,'1','N-D');
u(1)
trapz(x,u)
h = [0.04 0.02 0.01 0.005];
Eh = [];
u_ex = @(x) 4*cos(pi/2*x)+1-x;
for i = h
    [x,u]= miste(a,b,alpha,gamma,mu,0,sigma,f,(b-a)/i-1,'1','N-D');
    Eh = [Eh max(abs(u'-u_ex(x)))];
end
Eh
stimap_2(Eh,h);
%%
f = @(x) cos(3*x + sqrt(5));
a = 0; b = 1; poly_stima(f,4,a,b,'equi')
%%
f = @(x) exp(x); interp_tratti(f,-2,2,2,4,'equi',1.5)
%%
x = 0:10;
f = @(x) exp(x./10)+0.1*sin(pi*x+sqrt(2));
polyval(polyfit(x,f(x),2),11)
%%
gausslegendre_comp(0,4,@(x) exp(x),1,2,'equi')
%%
syms g;
f = @(x) g*x^3-4;
h = 1/3;
e = -1/12*h^2*12*g
%%
f = @(t,y) 3*exp(y)-19*t;
y0 = 0;
tf = 10;
tv = [0 tf]; h = 0.1;
[t,u] = eulero_avanti(f,10,y0,h);
u(2)
%%
mu = 1; sigma = @(x) 0.*x + 5; f = @(x) 0.*x +1; eta = 0;
a = 0; b = 1; alpha = 1; beta = 0; h = 0.1; N = (b-a)/h-1;
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
u(6)
%%
e1 = 8*1e-3; r = 2;
e2 = e1*(1/r)^2
%%
mu = 1; a = 0; b = 1;
u_s = @(t) 0.*t;
u_d = @(t) 0.*t;
g_0 = @(x) 3*sin(pi*x);
T = 1; delta_t = 1/8; h = 1/2; theta = 1;
f = @(t,x) 0.*x;
[u, x, t] = EqCalore_DiffFin_Theta(mu, f, a, b, u_s, u_d, g_0, T, h, delta_t, theta);
u(2,2)
%%
f = @(t) 10;
A = [-2 -6; 1 0];
g = @(t) [f(t) 0]';
tf = 5; tv = [0 tf]; y0 = [1 4]';
fun = @(t,y) A*y + g(t);
%
h = 0.1;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(2,2),u(2,end)
x = @(t) exp(-t).*(7/3.*cos(sqrt(5).*t)+10./(3.*sqrt(5)).*sin(sqrt(5).*t))+5/3;
h = [10 5 2.5 1.25].*1e-3;
Eh = [];
plot(t,u(2,:),t,x(t));
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh abs(u(2,end)-x(tf))];
end
Eh, stimap_2(Eh,h);
h = 0.1;
[t,u] = Heun_sistemi(fun,tv,y0,tf/h);
u(2,2),u(2,end)
l = max(eig(A));
fsolve(@(h) abs(1+h*l+(h*l)^2/2)-1,1)
%%
mu = 1; sigma = @(x) (1+19*x); f = @(x) 0.*x + 5; eta = 0;
a = 0; b = 1; alpha = 0; gamma = 0;
[x,u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,9,'1','N-N',[],0.1);
u(6)
%%
e1 = 2*1e-2; h1 = 0.1; h2 = 0.05; e2 = e1*(h2/h1) 
%%
A = [-4 1; 2 -3]; g = @(t) 0;
y0 = [7 2]';
[t,u]=eulero_indietro_sistemi_lineari(A,g,[0 1],y0,1/0.1);
u(:,end)
%%
f  = @(x) exp(3.*x); a = -1; b = 1;
poly_stima(f,3,a,b,'CGL')
%%
A = @(y) [-2 -10*y(2);1 0];
z = @(t) exp(-t/2)*(2*cos(t)-7/2*sin(t))+40*exp(-t)*sin(t)^2;
g = @(t) [z(t) 0]';

fun = @(t,y) A(y)*y + g(t); y0 = [2 0]';
tf = 10; tv = [0 tf]; 

h = 1e-2;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(2,2),u(2,end)
%%
A = diagonals([3 -2 -1],9);
g = @(t) exp(sin(pi*t))*ones(9,1);
fun = @(t,y) A*y + g(t);
y0 = 4*ones(9,1); tf = 10; tv = [0 tf];
h = 0.1;
[t,u1] = eulero_avanti_sistemi(fun,tv,y0,tf/h);

[t,u2] = Crank_Nicolson_sistemi_lineari(A,g,tv,y0,tf/h);
subplot(2,1,1)
plot(t,u1,'LineWidth',2)
subplot(2,1,2)
plot(t,u2,'LineWidth',2)
%%
clear
clc

z = @(t) exp(-t/2)*(2*cos(t)-7/2*sin(t))+40*exp(-t)*sin(t)^2;
f = @(y) [y(2); -2*y(2)-10*y(1)^2];
g = @(t) [0; z(t)];
y0 = [0 2]';
fun  = @(t,y) f(y)+g(t);
tf = 10; tv = [0 tf]; h = 1e-2;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
