%%
f = @(x) sin(x+sqrt(7));
a = 0; b = pi; n = 3;
[y_pol, P, err] = poly_lagr(f,n,a,b,pi,'equi');
dP = polyder(P); ddP = polyder(dP);
polyval(dP,pi)
polyval(ddP,pi)
%%
f = @(x) exp(2*(1-x));
a = 0; b = 2; n = 6;
poly_stima(f,n,a,b,'equi')
%%
f = @(x) x.^3;
interp_tratti(f,-1,1,2,4,'equi',0.75)
%%
x = 0:5;
y = [4 4 1 1 0 0];
poly_minim_IMQ(2,x,y)
%%
f = @(x) exp(x); a = -3; b = 1;
trapcomp(a,b,(b-a)/1,f)
%%
M1 = 25; e1 = 0.08; e2 = 0.02;
M2 = ceil(M1*nthroot(e1/e2,2))
%%
f = @(t,y) -2*y; y0 = 10; teta = 3/4; h = 0.1; tf = 10; tv = [0 tf];
[t,u] = teta_met(f,tv,y0,tf/h,teta);
u(2)
%%
mu = 1; eta = 0; sigma = @(x) 0.*x+3; f = @(x) -3.*x.^2.*(4-x.^2);
a = 0; b = 1; alpha = 0; beta = 1;
h = 0.1; N = 9;
u_ex = @(x) x.^4;
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
max(abs(u'-u_ex(x)))
%%
mu = 1; eta = 100; sigma = @(x) 0.*x; f = @(x) 0.*x;
a = 0; b = 1; alpha = 0; beta = 5; h = 0.1; N = 9;
Peh = abs(eta)*h/(2*mu);
mu = mu*(1+Peh);
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
u(10)
%%
mu = 1; f = @(x,t) 0.*x; a = 0; b = 1; u_s = @(t) 0; u_d = @(t) 0;
g_0 = @(x) 10*sin(pi*x); T = 2; delta_t = 0.1; h = 0.5; theta = 0;
[u, x, t] = EqCalore_DiffFin_Theta(mu, f, a, b, u_s, u_d, g_0, T, h, delta_t, theta);
u(2,4)
%%
A = diagonals([4 -2 -2],10);
g = @(t) cos(pi*t)*ones(10,1);
y0 = 7*ones(10,1);
fun = @(t,y) A*y + g(t);
tf = 10; tv = [0 tf]; h = 0.1;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
[u(3,2) u(3,end)]
[u_min,i] = min(u(3,:));
[u_min t(i)]

y_es = 0.269540864495346;
h = [0.05 0.025 0.0125 0.00625];
Eh = [];
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh abs(u(3,end)-y_es)];
end
p = stimap_2(Eh,h);
Eh
p(end)
loglog(h,Eh,h,h,h,h.^2,'LineWidth',2)
legend('Eh','h','h^2');
l = max(eig(A));
fsolve(@(h) abs(1+h*l)-1,1)
h = 0:0.001:1;
plot(h,abs(1+h.*l),'LineWidth',2);
yline(1,'LineWidth',2);

[t,u] = eulero_indietro_sistemi_lineari(A,g,tv,y0,tf/0.1);
[u(3,2) u(3,3)]
A = [1/4 0; 1/2 1/4];
b = [1/2 1/2];
c = [1/4 3/4];
[t,u] = Runge_Kutta(A,b,c,fun,tv,y0,tf/0.1);
[u(3,2) u(3,end)]
%%
x = 0:0.5:1; y = [2 1 1.5];
[P,S] = polyfit(x,y,2);
[y,delta]=polyval(P,0.25,S); y
polyval(polyder(P),0.25)
%%
x = 0:0.25:1; y = [3 0.5 1.5 -0.5 1];
x_dis = 0:0.01:1;
min(interp1(x,y,x_dis))
%%
f = @(x) sin(pi.*x);
n = 4;
h = 2/n;
x_nod = 0:h:2;
interp1(x_nod,f(x_nod),1.6)
%%
poly_scarto_IMQ(0:0.25:1,[3 0.5 1.5 -0.5 1],2)
%%
f = @(x) sqrt(2+abs(x));
N = 3;
trapcomp(-1,2,N,f)
%%
p = 2;
M1 = 20; e1 = 1e-1; e2 = 1e-3;
M2 = ceil(M1*nthroot(e1/e2,p))
%%
f = @(t,y) -(1+sin(t))*y.^2./81; h = 0.2; y0 = 9;
df = @(t,y) -2*(1+sin(t))*y/81;
[t,u] = eulero_indietro_newton(f,df,5,y0,h);
u(end)
%%
mu = 1; eta = 0; sigma = @(x) 2 + 0.*x; f = @(x) 3 + x;
a = 0; b = 1; alpha = 1; beta = 0; h = 0.1; N = 9;
[x,u] = miste(a,b,alpha,beta,mu,eta,sigma,f,N,'1','D-N','',h);
u(end)
%%
a = 0; b = 1; alpha = 0; beta = 3; mu = 1; eta = -50; sigma = @(x) 0.*x; f = @(x) 0.*x;h = 0.1; N = 9;
Peh = abs(eta)*h/(2*mu);
mu = mu*(1+Peh);
[x_nodes, u] = Dirichlet(a, b, alpha, beta, mu, eta, sigma, f, N, h);
u(2)
%%
K = 1; h = 1e-2; 
A = @(y) [-1 -(147/16+K*y(2)); 1 0];
g = @(t) [-9*exp(-t/2)*(sin(3*t)^2-1)-9/2*exp(-t/4)*sin(3*t);0];
fun = @(t,y) A(y)*y + g(t);
tf = 5; y0 = [-3/4 3]'; tv = [0 tf];

[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(:,2),u(:,end)
[u2_min,i] = min(u(2,:));
u2_min, t(i)


y = @(t) 3*exp(-t/4).*[-1/4*cos(3*t)-3*sin(3*t); cos(3*t)];
h = [10 5 2.5 1.25].*1e-1; Eh = [];
for i = h
    [t,u] = Crank_Nicolson_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh norm(u(:,end)-y(t(end)))];
end
Eh
stimap_2(Eh,h);
loglog(h,Eh,'-',h,h.^2,'-.',h,h,':','LineWidth',2)
legend('Eh','h^2','h')
%%
RK = [0 0;2/3 0]; b = [1/4 3/4]; c = [0 2/3]; h = 1e-2;
[t,u] = Runge_Kutta_gen(RK,b,c,fun,tf,h,y0);
u(:,end)
%%
K = 0; h = 1e-2; 
A = [-1 -(147/16); 1 0];
fun = @(t,y) A*y + g(t);
tf = 5; y0 = [-3/4 3]'; tv = [0 tf];
g = @(t) zeros(2,1).*t;
[t,u] = Crank_Nicolson_sistemi_lineari(A,g,tv,y0,tf/h);
u(:,2),u(:,end)
%%
f = @(x) sin(x+sqrt(5));
x_dis = 0:0.001:pi; n = 30; a = 0; b = pi;
k = 0:n;
    t = -cos(pi*k/n);
    x_nod = ((b-a)/2)*t + (a+b)/2;

[y,P] = poly_lagr(f,n,0,pi,x_dis,'CGL');

dP = polyder(P);
y2 = polyval(dP,x_dis);
plot(x_nod,f(x_nod),'o',x_dis,y,x_dis,y2)
%%
mu = 2; eta = 0; sigma = @(x) 1 + 0.*x; f = @(x) 3*sin(x);
a = 0; b = 1; alpha = 0; beta = sin(1); N = 9; h = 0.1;
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
u_ex = @(x) sin(x);
E = max(abs(u'-u_ex(x)))
%%
A = [-1/6 -5; 1 0]; 
z = @(t) 0;
g = @(t) [z(t)/2 0]';
tf = 50; tv = [0 tf];
fun = @(t,y) A*y + g(t);
y0 = [-3/4 9]'; h = 0.1;
[t,u] = Heun_sistemi(fun,tv,y0,tf/h);
u(2,2),u(2,end)
%%
m = 10;
A = diagonals([4 -2 -2],m);
g = @(t) cos(pi*t).*ones(m,1);
fun = @(t,y) A*y + g(t);
y0 = 7*ones(m,1);
tf = 10; tv = [0 tf]; h = 0.1;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(3,2),u(3,end)
[t,u] = eulero_indietro_sistemi_lineari(A,g,tv,y0,tf/h);
u(3,2),u(3,3)

RK = [1/4 0; 1/2 1/4];
b = [1/2 1/2]; c = [1/4 3/4];
[t,u] = Runge_Kutta_linear(RK,b,c,A,g, tv, y0, tf/h);
u(3,2),u(3,end)
%%
clear
clc

A = [-4 -100; 1 0]; 
g = @(t) [200*cos(10.*t) 0]';
y0 = [50 0]'; tf = 5; tv = [0 tf];
fun = @(t,y) A*y + g(t);

RK = [0 0 0; 1/2 0 0; -1 2 0];
b = [1/6 2/3 1/6]; c = [0 1/2 1]; h = 1e-2;
[t,u] = Runge_Kutta_linear_espl(RK,b,c,A,g,tv, y0, tf/h);
u(2,end)
plot(t,u)
%%
clear
clc
m = 9;
A = diagonals([5/2 -2 -1/2],m);
g = @(t) (1+2*sin(pi*t))*ones(m,1);
y0 = 2*ones(m,1);
tf = 10; tv = [0 tf]; h = 0.1;
RK = [1/4 0; 1/2 1/4];
b = [1/2 1/2];
c = [1/4 3/4];
[t,u] = Runge_Kutta_linear_impl(RK,b,c,A,g,tv,y0,tf/h);
u(5,2),u(5,end)
%%
clear
clc
a = -1; b = 1;
x = a:0.001:b; 
n = 5;
i = 0:n;
x_nod = -cos(pi.*i./n);
syms y;
omega = 1;
for i = 1 : n + 1
omega = omega*(y-x_nod(i));
end
omega = matlabFunction(simplify(expand(omega)));
stim = (1e3)*max(abs(omega(x)))/factorial(n+1)
%%
clear
clc
f = @(t,y) -(t+1)*y^2; t0 = 0; y0 = 9;
tf = 10; tv = [t0 tf];
RK = [0 0; 1 0]; b = [1/2 1/2]; c = [0 1]; h = 0.1;
[t,u]  = Runge_Kutta(RK,b,c,f,tv,y0,tf/h);
u(2)
%%
syms t y;
f =  -t*y;
HN_sym(f,1,5)
%%
clear
clc

V = 1; h = 0.1;
mu = 1; eta = V; sigma = @(x) 0.*x; f = @(x) sin(x)+cos(x);
a = 0; b = 1; alpha = 0; beta = cos(1);
N = (b-a)/h - 1;
[x,u] = miste(a,b,alpha,beta,mu,eta,sigma,f,N,'1','D-N','',h);
u(end)
I = 0;
for i = 2 : length(u)
     I = I + abs((u(i-1)+u(i))*h/2);
end
I 


h = [0.1 0.05 0.025 0.00125];
u_ex = @(x) sin(x); Eh = [];
for i = h
    N = (b-a)/i - 1;
    [x,u] = miste(a,b,alpha,beta,mu,eta,sigma,f,N,'1','D-N','',i);
    Eh = [Eh max(abs(u'-u_ex(x)))];
end
stimap_2(Eh,h);

clear
clc

V = 100; h = 0.1;
mu = 1; eta = V; sigma = @(x) 0.*x; f = @(x) 0.*x;
a = 0; b = 1; alpha = 0; beta = cos(1);
N = (b-a)/h - 1;
Peh = abs(eta)*h/(2*mu)

[x,u] = miste(a,b,alpha,beta,mu,eta,sigma,f,N,'1','D-N','up',h);
plot(x,u)
u(end)
%%
clear
clc
m = 9;
A = diagonals([5/2 -2 -1/2],m);
g = @(t) (1+2*sin(pi*t))*ones(m,1);
y0 = 2*ones(m,1);
tf = 10; tv = [0 tf]; 

h = [10 5 2.5 1.25].*1e-2;
Eh = [];
for i = h
    [t,u] = Crank_Nicolson_sistemi_lineari(A,g,tv,y0,tf/i);

end
%%
n = 4; a = 0; b = 4; h = (b-a)/n;
x_nod = a:h:b;
polyval(polyfit(x_nod,[-4 -2 -2 0 -2],n),2.5)
%%
f = @(x) 3-7*x^3;
pmedcomp(-4,4,1,f)
%%
f = @(x) exp(x);
integr_stima(0,3,1e-3,f,'t')
%%
f = @(x) 1-3^(8*x);
h = 1/8;
x = 0;
df = (f(x)-f(x-h))/h
f = @(x) 3-5*x^2;
%%
mu = 1; f = @(x) 6 + 0.*x; sigma = @(x) 0.*x; eta = 0;
a = 0; b = 1; alpha = 0; beta = 0;
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,1,0.5);
u(2)
%%
mu = 1; eta = 0; sigma = @(x) 0.*x + 5;
f = @(x) 8*sin(pi.*x);
a = 0; b = 1; alpha = 0; N = 1; beta = 1; h = 0.5;
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
u(2)
%%
f = @(x) 8.*(exp(x./10)+sin(pi*x+sqrt(3)));
a = 0; b = 10;
H = [0.25 0.125 0.0625 0.03125];
Eh = [];
for i = H
    x_nod = a:i:b;
    interp1(x_nod,f(x_nod),sqrt(2))
    x_dis = a:0.001:b;
    Eh = [Eh max(abs(interp1(x_nod,f(x_nod),x_dis)-f(x_dis)))];
end
Eh
n = 10;
h = (b-a)/n;
x_nod = a:h:b;
polyfit(x_nod,f(x_nod),1)
polyfit(x_nod,f(x_nod),2)
%%
f = @(t,y) 5*(3*cos(3*t)*exp(-t/2)+1/2)-y/2;
h = [0.5 0.25 0.125 0.0625]; tf = 3; y0 = 5;
U = []; Eh = [];
y = @(t) 5*(1+sin(3.*t).*exp(-t./2));
for i = h
    [t,u] = Crank_Nicolson(f,tf,y0,i);
    U = [U u(end)];
    Eh = [Eh max(abs(u(end)-y(tf)))];
end
U 
Eh
%%
a = 0; b = 1; c = 0; d = 1;
x = [0.25 0.75];
y = x;
S = 0;
M = 4;
f = @(x,y) exp(2*x+y);
for i = 1:2
    for j = 1:2
        S = S + f(x(i),y(j));
    end
end
I = (b-a)*(d-c)*S/M
%%
a = 0; b = pi;
f = @(x) 5.*(x.^3+sin(4*x+sqrt(2)));
pmedcomp(a,b,1,f)
pmedcomp(a,b,10,f)

%%
n = 2;
y = [-sqrt(3/5) 0 sqrt(3/5)];
alpha = [5/9 8/9 5/9];
gausslegendre_comp(0,pi,f,1,2,'equi')
gausscomp_2(a,b,1,f)
%%
f = @(x) 5*(1-x.^4+3*x.^2);
a = 0; b = 1;
gausslegendre_comp(a,b,f,1,1,'equi')
%%
f = @(x) 3*x.^2.*(1+x.^2);
gausslegendre_comp(0,1,f,1,1,'equi') 
%%
a = 0; b = 10; f = @(x) (x+2).^2;
n = 5;  h = (b-a)/n;
x_nod = a:h:b;
interp1(x_nod,f(x_nod),1)
%%
f = @(x) sin(10.*x)-7.*x;
a = 0; b = 3; tol = 1e-4; 
num_int = interpH(f,tol,a,b)
%%
clear
clc
V = 1; h = 0.1; a = 0; b = 1; alpha = 1; gamma = exp(1);
mu = 1; eta = V; sigma = @(x) 1 + 0.*x; f = @(x) exp(x);
N = (b-a)/h-1;
[x,u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,'1','D-N','',h);
u(end)
I = 0;
for i = 2:length(u)
    I = I+(u(i-1)+u(i))*h/2;
end
I
trapz(x,u)
%%
f = @(x) exp(x);
a = -2; b = 2;
tol = 1e-1;
integr_stima(a,b,tol,f,'t')
%%
f = @(t,y) -2*(t+1)*y^2; tf = 10;
tv = [0 tf]; y0 = 3;
A = [0 0; 3/4 0]; b = [1/3 2/3]; c = [0 3/4];   
[t,u] = Runge_Kutta(A,b,c,f,tv,y0,tf/h);
format long
u(2)
format short
%%
f = @(t,y) -5*(1+t)*y; tf = 100;
tv = [0 tf]; y0 = 8; h = 0.2;
[t,u] = Crank_Nicolson(f,tf,y0,h);
u(2)
%%
A = [-3 -4; 1 0]; l = max(eig(A))
fsolve(@(h) abs(1+h.*l)-1,1)
%%
mu = 3; eta = -50;
2*mu/abs(eta)
%%
format long
n = 5; a = 0; b = 1;
h = (b-a)/n;
x_nod = a:h:b;
%%
format long
f = @(x) exp(x);
a = -2; b = 2; H = 1; 
[y] = interp_tratti(f,a,b,2,4,'equi',1.5)
%%
x = 0:10;
f = @(x) exp(x./10) + 0.1*sin(pi.*x+sqrt(2));
polyval(polyfit(x,f(x),2),11)
%%
f = @(x) exp(x);
a = 0; b = 4;
gausslegendre_comp(a,b,f,1,2,'equi')
%%
clear
clc
mu = 1; eta = 0; sigma =@(x) 0.*x + 5; f = @(x) 0.*x +1;
a = 0; b = 1; alpha = 1; beta = 0;
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,9,0.1);
u(6)
%%
e1 = 8*1e-3; h1 = 1; h2 = 0.5;
e2 = e1*(h2/h1)^2
%%
format short
clear
clc
A = [-2 -6; 1 0];
y0 = [1 4]';
g = @(t) [10; 0];
fun = @(t,y)A*y + g(t);
tf = 5; tv = [0 tf]; h = 0.1;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(2,2),u(2,end)
x_es = @(t) exp(-t)*(7/3*cos(sqrt(5)*t)+10/(3*sqrt(5))*sin(sqrt(5)*t))+5/3;
h = [10 5 2.5 1.25].*1e-3; Eh = []; 
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh abs(u(2,end)-x_es(tf))];
end
Eh, stimap_2(Eh,h);
[t,u] = Heun_sistemi(fun,tv,y0,tf/0.1);
u(2,2),u(2,end)
l = max(eig(A));
fsolve(@(h) abs(1+h*l+(h*l)^2/2)-1,1)
h = 0.5
abs(1+h*l+(h*l)^2/2)
%%
hl = -4:0.001:4;
R = 1 + hl + hl.^2./2 + hl.^3/6 + hl.^4./factorial(4);
plot(hl,abs(R));yline(1)
% -2.7853 < hl < 0
% l < 0 --> 0 < h < -2.7853 / |l|
%%
clear 
clc
f = @(x) exp(2*x+sin(pi*x)); a = -1; b = 1; n = 5;
h = (b-a)/n; x_nod = a:h:b; x_dis = a:0.001:b;
y_dis = polyval(polyfit(x_nod,f(x_nod),n),x_dis);
[E,i] = max(abs(y_dis-f(x_dis)));
E, x_dis(i)
%%
x = 0:0.25:1; y = [1 0.5 1.5 -0.25 1];
plot(x,y,'o')
%%
f = @(x) 2*abs(sin(pi*x));
interp_tratti(f,0,4,2,4,'equi',1.75)
%%
phi = poly_phi([0 0.5 2]);
phi = matlabFunction(phi(1));
simpcomp(0,2,1e5,phi)
%%
M1 = 10; e1 = 1e-1; M2 = 100; e2 = 1e-4;
p = log(e1/e2)/log(M2/M1)
%%
a = 0; b = 1; c = 0; d = 1;
f = @(x,y) 2.^(x+3.*y);
I = (b-a)*(d-c)*(f(a,c)+f(a,d)+f(b,c)+f(b,d))/4
integral2(f,a,b,c,d)
%%
f = @(t,y) -sqrt(10*y/(10+t));
tf = 4; tv = [0 4]; y0 = 9;
h = 0.2; 
[t,u] = Heun(f,tf,y0,h);
u(end)
%%
h1 = 0.1; e1 = 2*1e-2; h2 = 0.05;
e2 = e1*(h2/h1)
%%
clear
clc
mu = 1; eta = 0; sigma = @(x) 0.*x + 3; f = @(x) 10.*sin(pi.*x);
a = 0; b = 1; alpha = 1; beta = 0;
h = 0.1;N = 9;
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
u(3)
%%
mu = 1; f = @(t,x) 0; a = 0; b = 1; u_s = @(t) 0; u_d = @(t) 0; g_0 = @(x) 6*sin(pi*x); T = 4; h = 0.5; delta_t = 0.1;
theta = 0.5;
[u, x, t] = EqCalore_DiffFin_Theta(mu, f, a, b, u_s, u_d, g_0, T, h, delta_t,theta);
%%
A = [-4 -100; 1 0]; g = @(t) [200*cos(10*t);0];
fun = @(t,y) A*y + g(t);
tf = 5; tv = [0 tf]; y0 = [50 0]';
h = 1e-2;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(2,2),u(2,end)
[u_min,i] = min(u(2,:)); u_min,t(i)
x_es = @(t) 5*sin(10*t); h = [10 5 2.5 1.25].*1e-4; Eh = [];
x_fin = x_es(tf);
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh abs(u(2,end)-x_fin)];
end
Eh,stimap_2(Eh,h);
%%
l = max(eig(A));
h = 0:0.01:1;
R = 1 + h.*l;
plot(h,abs(R),'LineWidth',2),yline(1,'LineWidth',2)
options = optimoptions('fsolve','Display','iter');
options = optimoptions('fsolve','Display','off');
fsolve(@(h) abs(1+h*l)-1,1,options)
%%
clear
clc
tf = 5; tv = [0 tf]; y0 = [50 0]';
h = 1e-2;
A = [-4 -100; 1 0]; g = @(t) [200*cos(10*t);0];
RK = [0 0 0; 0.5 0 0; -1 2 0];
b = [1/6 2/3 1/6]; c = [0 0.5 1];
[t,u] = Runge_Kutta_linear_espl(RK,b,c,A,g,tv,y0,tf/h);
u(2,2),u(2,end)
%%
f = @(x) exp(2*x); n = 5; a = -1; b = 1; x_dis = a:0.0001:b;
[y_pol, P, err] = poly_lagr(f,n,a,b,x_dis,'equi');
err
poly_stima(f,n,a,b,'equi')
%%
x = 1:8;
y = [679 776 882 794 932 808 480 907]; x_dis = 1:0.01:12;
y_dis = polyval(polyfit(x,y,2),x_dis);
plot(x,y,'o',x_dis,y_dis,'LineWidth',2)
%%
f = @(x) exp(x/2);
interp_tratti(f,0,3,2,3,'equi',0.75)
%%
f = @(x) exp(x); a = -1; b = 2; trapcomp(a,b,1,f)
%%
M1 = 10; e1 = 1e-1; M2 = 100;
p = 2; e2 = e1*(M1/M2)^p
%%
syms t h y;
f = -sqrt(17*y/(17+t^2));
EA_sym(f,1,9)
%%
h1 = 1e-3; e1 = 4*1e-3; h2 = 5*1e-4;
e2 = e1*(h2/h1)^2
%%
mu = 1; eta = 40; sigma = @(x) 0.*x; f = @(x) 0.*x;
a = 0; b = 1; alpha = 3; beta = 0; h = 0.1; N = 9;
Peh = abs(eta)*h/(2*mu); mu = mu*(1+Peh);
[x_nodes,uh] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
uh(10)
%%
clear
clc
m = 9;
A = diagonals([5/2 -2 -0.5],m); g = zeros(m,1);

l = max(eig(A));
options = optimoptions('fsolve','Display','off'); 
h_max = fsolve(@(h) abs(1+h*l)-1,1,options)

tf = 10; g = @(t) (1+2*sin(pi*t))*ones(m,1); y0 = 2*ones(m,1);
h = 0.1;
fun = @(t,y) A*y + g(t); tv = [0 tf];
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(5,2),u(5,end)
[u_min,i] = min(u(5,:)); u_min,t(i)
y_es = 1.01702435;
h = [10 5 2.5 1.25].*1e-3;
Eh = [];
for i = h 
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);  
    Eh = [Eh abs(u(5,end)-y_es)];
end
Eh, stimap_2(Eh,h); h = 0.1;
[t,u] = eulero_indietro_sistemi_lineari(A,g,tv,y0,tf/h);
u(5,2)
RK = [1/4 0; 1/2 1/4]; b = [0.5 0.5]; c = [1/4 3/4];
[t,u] = Runge_Kutta_linear_impl(RK,b,c,A,g,tv,y0,tf/h);
u(5,2),u(5,end)
%%
f = @(x) sqrt(x+sin(2*pi*x)+3); a = -1; b = 1; 
poly_lagr(f,5,a,b,0.4,'equi')
%%
n = 5; a = -1; b = 1;
der = 1e3;
x = a:0.001:b;
i = 0:n;
x_nod = -cos(pi.*i./n);
syms y;
omega = 1;
for i = 1 : n + 1
    omega = omega*(y-x_nod(i));
end
omega = matlabFunction(simplify(expand(omega)));
stim = max(abs(der))*max(abs(omega(x)))/factorial(n+1)
%%
f = @(x) exp(x/6); a = 0; b = 3; interp_tratti(f,a,b,4,3,'equi',1.5)
%%
f = @(x) 1 + x.^5 + abs(x);
trapcomp(-1,1,1,f)
%%
M1 = 20; e1 = 0.1; M2 = 40;
e2 = e1*(M1/M2)^2
%%
Pi = [0 0; 1 0;0 1];
f = @(x) exp(2*x(1)+x(2));
T = 0.5;
P = sum(Pi,1)/3;
I = T*f(P)
%%
syms t y;
f = -t*y; y0 = 5;
EI_sym(f,1,y0)
%%
A = [0 0; 1 0]; b = [0.5 0.5]; c = [0 1];
h = 0.1;
f = @(t,y) -(t+1)*y^2; t0 = 0; tf = 10; y0 = 9; tv = [0 tf];
[t,u] = Runge_Kutta(A,b,c,f,tv,y0,tf/h);
u(2)
%%
A = [-2 1; 2 -3]; y0 = [9 3]';
h = 0.1; fun = @(t,y) A*y;
[t,u] = eulero_indietro_sistemi(fun,[0 1],y0,1/h);
u(:,end)
%%
f = @(t,u,v) -2*u; y0 = 10; w0 = 0; t_vett = [0 10]; h = 0.1;
[t, u, v] = Leap_Frog(f, t_vett, h, y0, w0);
u(2),u(3)
%%
V = 1;
mu = 1; eta = V; sigma = @(x) 0.*x; f = @(x) sin(x)+cos(x);
a = 0; b = 1; alpha = 0; gamma = cos(1);
h = 0.1; N = (b-a)/h-1;
[x,u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,'1','D-N','',h);
u(end)
%%
I = 0;
for i = 2: length(u)
    I = I + (u(i-1)+u(i))*h/2;
end
I
h = [0.1 0.05 0.025 0.00125]; Eh = []; u_ex = @(x) sin(x);
for i = h
    N = (b-a)/i-1;
    [x,u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,'1','D-N','');
    Eh = [Eh max(abs(u'-u_ex(x)))];
end
Eh, stimap_2(Eh,h);
V = 100;53-15
mu = 1; eta = V; sigma = @(x) 0.*x; f = @(x) 0.*x;
a = 0; b = 1; alpha = 0; gamma = cos(1);
h = 0.1; N = (b-a)/h-1;
[x,u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,'1','D-N','up',h);
u(end)
%%
f = @(x) exp(x).*(1+sin(pi.*x));
n = 5; a = -1; b = 1; x_dis = a:0.001:b;
[y_pol, P, err] = poly_lagr(f,n,a,b,x_dis,'equi');
err
%%
T = [6 5.5 6 4.5 4 3.5 4 4.5 5];
g = 1:9;
polyval(polyfit(g,T,2),10)
%%
syms x g;
f = x^g;
df = jacobian(f,x)
ddf = jacobian(df,x)
%%
f = @(x) sqrt(1+x); a = -1; b = 3; trapcomp(a,b,4,f)
%%
M1 = 10; e1 = 1e-1; M2 = 100; e2 = e1*(M1/M2)^4
%%
syms t h y;
f = -(1+t)*y; y0 = 2;
HN_sym(f,1,y0)
%%
mu = 1; sigma = @(x) 0.*x; f = @(x) (2+x).^2; eta = 0;
a = 0; b = 1; alpha = 1; beta = 0; N = 9; h = 0.1;
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
u(6)
%%
e1 = 4*1e-3; h1 = 1e-4; h2 = 5*1e-3;
e2 = e1*(h2/h1)
%%
mu = 1; eta = 40; sigma = @(x) 0.*x; f = @(x) 0.*x;
a = 0; b = 1; alpha = 7; gamma = 0; h = 0.1; N = 9;
Peh = abs(eta)*h/(2*mu);
mu = mu*(1+Peh);
[x,u] = Dirichlet(a,b,alpha,gamma,mu,eta,sigma,f,N,h);
u(10)
%%
x = 0:4; y = [3 1 3 4 5];
polyval(polyfit(x,y,4),1.5)
%%
f = @(x) exp(x);
poly_stima(f,4,0,2,'equi')
%%
f = @(x) x.^3+sin(pi.*x);
interp_tratti(f,0,3,2,3,'equi',1.5)
%%
f = @(x) x.^2+abs(x);
a = -1; b = 1; x_dis = a:0.0001:b;
n = 4;
h = (b-a)/n; x_nod = a:h:b;
y_dis = polyval(polyfit(x_nod,f(x_nod),2),x_dis);
trapz(x_dis,y_dis)
%%
M1 = 10; e1 = 1e-1; M2 = 100;
e2 = e1*(M1/M2)^2
%%
syms t h y;
f = -t^2*y; y0 = 1;
EI_sym(f,1,y0)
%%
A = [0 0; 0.5 0]; b = [0 1]; c = [0 0.5];
f = @(t,y) -(t+1)*y^3; t0 = 0; tf = 10; tv = [0 tf]; y0 = 1;
h = 0.1;
[t,u] = Runge_Kutta(A,b,c,f,tv,y0,tf/h);
u(2)
%%
A = [-2 1; 2 -3]; y0 = [1 3]'; tf = 1; tv = [0 tf];
fun = @(t,y) A*y;
[t,u] = Crank_Nicolson_sistemi(fun,tv,y0,tf/h);
u(:,end)
%%
A = [-2 -1; 1 0]; y0 = [1 4]';
tf = 10; tv = [0 tf]; h = 0.1;
fun = @(t,y) A*y;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(2,2)
%%
f = @(t,u,v) -u;
y0 = 5; w0 = 0; h = 0.1; 
[t,u] = Leap_Frog(f,[0 10],h,y0,w0);
u(2:4)
%%
V = 1;
mu = 1; eta = V; sigma = @(x) 0.*x + 1; f = @(x) exp(x);
a = 0; b = 1; alpha = 1; gamma = exp(1);
h = 0.1; N = (b-a)/h-1;
[x,u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,'1','D-N','');
u(2),u(end)
I = 0;
for i = 2:length(u)
    I = I + (u(i-1)+u(i))*h/2;
end
I
h = [0.1 0.05 0.025 0.0125];
u_ex = @(x) exp(x); Eh = [];
for i = h
    N = (b-a)/i-1;
    [x,u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,'1','D-N','');
    Eh =[Eh max(abs(u'-u_ex(x)))];
end
Eh, stimap_2(Eh,h);
V = 1000;
mu = 1; eta = V; sigma = @(x) 0.*x + 1; f = @(x) 0.*x;
a = 0; b = 1; alpha = 1; gamma = exp(1);
h = 0.1; N = (b-a)/h-1;
[x,u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,'1','D-N','up');
u(end)
%%
f = @(x) exp(x).*(1+sin(pi*x));
x = -1:0.0001:1;
[y_pol, P, err] = poly_lagr(f,5,-1,1,x,'equi');
err
%%
g = 1:9;
T = [6 5.5 6 4.5 4 3.5 4 4.5 5];
polyval(polyfit(g,T,2),10)
%%
syms x g;
f = x^g;
df = jacobian(f,x);
ddf = jacobian(df,x);
h = 1/10;
e = h^2/8*subs(ddf,x,1)
%%
f = @(x) sqrt(1+x);
trapcomp(-1,3,4,f)
%%
M1 = 10; e1 = 1e-1; M2 = 100;
e2 = e1*(M1/M2)^4
%%
clear
clc
a = 0; c = 0; b = 1; d = 1;
eps = [-1/sqrt(3) 1/sqrt(3)];
x = (a+b)/2+(b-a)/2.*eps;
y = (c+d)/2+(d-c)/2.*eps;
I = 0;
f = @(x,y) exp(x+3*y);
for i = 1:2
    for j = 1:2
        I = I + f(x(i),y(j));
    end
end
I = I*(b-a)*(d-c)/4
%%
syms t h y;
f = -(1+t)*y; y0 = 2;
HN_sym(f,1,y0)
%%
mu = 1; eta = 0; sigma = @(x) 0.*x; f = @(x) (2+x).^2;
a = 0; b = 1; alpha = 1; beta = 0;
h = 0.1; N = 9;
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
u(6)
%%
h1 = 1e-4; e1 = 4*1e-3; h2 = 5*1e-3;
e2 = e1*(h2/h1)
%%
mu = 1; eta = 40; sigma = @(x) 0.*x; f = @(x) 0.*x;
h = 0.1; N = 9; a = 0; b = 1; alpha = 7; beta = 0;
Peh = abs(eta)*h/(2*mu); mu = mu*(1+Peh);
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
u(10)
%%
m = 9;
A = diagonals([3 -2 -1],m); 

% con g(t) = 0 :
l = max(eig(A));
h_max = -2*real(l)/abs(l)^2
options = optimoptions('fsolve','Display','off');
fsolve(@(h) abs(1+h*l)-1,1)

g = @(t) exp(sin(pi*t))*ones(m,1); y0 = 4*ones(m,1);
fun = @(t,y) A*y + g(t);
h = 0.1; tf = 10; tv = [0 tf];
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(5,2),u(5,end)
[u5_min,i] = min(u(5,:));
u5_min, t(i)
y_es = 1.142174435142178;
h = [10 5 2.5 1.25].*1e-3;
Eh = [];
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh abs(u(5,end)-y_es)];
end
Eh,stimap_2(Eh,h);
loglog(h,Eh,h,h,'-.',h,h.^2,'-',LineWidth=2)
legend('Eh','H','H^2')
h = 0.1;
[t,u] = Crank_Nicolson_sistemi_lineari(A,g,tv,y0,tf/h);
u(5,2)
RK = [1/3 0; 1/3 1/3]; b = [0.5 0.5]; c = [1/3 2/3];
[t,u] = Runge_Kutta_linear_impl(RK,b,c,A,g,tv,y0,tf/h);
u(5,2),u(5,end)
%%
f = @(x) (x/pi).*(x<=pi) + ((2-x/pi).^2).*(x>pi);
a= 0; b = 2*pi;
n = 4;
x_nod = a:(b-a)/n:b;
y_nod = f(x_nod);
Y = fft(y_nod);
C = fftshift(Y)/(n+1);
x_dis = linspace(a,b,100);
y_dis = interpft(y_nod,100);
trig = matlabFunction(poly_trigo(f,n));
plot(x_nod,f(x_nod),'o',x_dis,trig(x_dis),x_dis,y_dis)

%%
fun = @(t,y) [-2*y(1).^2-10*y(2); y(1)];
x = @(t) 2*exp(-t/2)*sin(t);
J = @(t) [-4*y(1) -10; 1 0];
t = 0:0.001:10;
l = [];
for i = 1 : length(t)
    l = [l max(eig(J(t(i))))];
end
l = max(l);
%% 
clear
clc
A = [0 0;3/4 0]; b = [1/3 2/3]; c = [0 3/4];
f = @(t,y) -2*(t+1)*y^2; t0 = 0; tf = 10; tv = [0 tf]; y0 = 3; h = 0.1;
[t,u] = Runge_Kutta(A,b,c,f,tv,y0,tf/h);
format long
u(2)
format short
%%
f = @(x) sin(x+sqrt(2)); a = 0; b = pi; n = 3; x = 2;
[y_pol, P, err] = poly_lagr(f,n,a,b,x,'equi');
y_pol, polyval(polyder(P),2)
%%
f = @(x) exp(3.*x); a = -1; b = 1; n = 3;
poly_stima(f,n,a,b,'CGL')
%%
n = 5;
f = @(x) ((x./pi).^2).*(x>=0 & x<pi)+((2-x./pi).^2).*(x>=pi & x<2*pi);
x = 2*pi/n.*(0:n);
cubicspline(x,f(x),x(4))
%%
x = 0:4; y = [2 2 0 1 0];
poly_minim_IMQ(1,x,y)
%%
f = @(x) exp(x); a = -1; b = 3; trapcomp(a,b,4,f)
%%
f = @(x) x.^3;
I = (6)^4/4;
fsolve(@(c) 3*(f(3+c)+f(3-c))-I,1)
%%
f = @(t,y) -2*y; teta = 3/4;
h = 0.1; y0 = 7; tf = 10; tv = [0 tf];
[t,u] = teta_met(f,tv,y0,tf/h,teta);
u(2)
%%
h1 = 1e-2; e1 = 4*1e-3; e2 = 1e-3;
h2 = h1*nthroot(e2/e1,2)
%%
mu = 1; eta = 0; sigma = @(x) 1+19*x; f = @(x) 0.*x + 5; a = 0; b = 1;
alpha = 0; beta = 0;
h = 0.1; N = 9;
[x,u] = miste(a,b,alpha,beta,mu,eta,sigma,f,N,'1','N-N','');
u(6)
%%
mu = 1; eta = 0; sigma = @(x) 0.*x; f = @(x) 0.*x + 6; 
a = 0; b = 1; b_b = 3; alpha = 0; beta = 0; N = 9; h = 0.1;
[t,u] = Robin(a,b,alpha,beta,b_b,mu,eta,sigma,f,N,h,'D-R');
u(end)
%%
f = @(y) [-2*y(1)-10*y(2)^2;y(1)];
z = @(t) exp(-t/2)*(2*cos(t)-7/2*sin(t))+40*exp(-t)*sin(t)^2;
g = @(t) [z(t);0];
fun = @(t,y) f(y)+g(t); y0 = [2 0]';
tf = 10; tv = [0 tf]; h = 1e-2;
[t,u]=eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(2,2),u(2,end)
plot(t,abs(u(2,:)),'-o'),yline(0.1)
t_min = 5.68;

x_es = @(t) 2.*exp(-t./2).*sin(t);
Eh = [];
h = [10 5 2.5 1.25].*1e-4;
for i = h
    [t,u]=eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh max(abs(u(2,:)-x_es(t)))];
end
Eh, stimap_2(Eh,h);
loglog(h,Eh,h,h,'-',h,h.^2,':','LineWidth',2);
legend('Eh','H','H^2');
J = @(t) [-2 -20*x_es(t); 1 0];
t = 0:0.01:tf;
l = [];
for i = 1:length(t)
    l = [l max(eig(J(t(i))))];
end
l = max(l);
fsolve(@(h) abs(1+h*l)-1,1)
[t,u] = multipasso(fun,tv,y0,tf/1e-2);
u(2,2),u(2,3),u(2,end)
h = [10 5 2.5 1.25].*1e-4; Eh = [];
for i = h
    [t,u]=multipasso(fun,tv,y0,tf/i);
    Eh = [Eh max(abs(u(2,:)-x_es(t)))];
end
Eh, stimap_2(Eh,h);
%%
clear
clc
V = 1;
mu = 1; eta = V; sigma = @(x) 0.*x; f = @(x) sin(x)+cos(x);
a = 0; b = 1; alpha = 0; beta = cos(1);
h = 0.1; N = (b-a)/h-1; 
[x,u] = miste(a,b,alpha,beta,mu,eta,sigma,f,N,'1','D-N','');
u(end)
I = 0;
for i = 2 : length(u)
    I = I + (u(i-1)+u(i))*h/2;
end
I 
%%
f = @(x) exp(x);
a = [5/9 8/9 5/9];
x = [-sqrt(3/5) 0 sqrt(3/5)];
I = sum(a.*f(x))
gausslegendre_comp(0,4,f,1,2,'equi')
%%
f = @(x) x.^3+3; a = 0; b = 10; n = 5;
x_nod = a:(b-a)/n:b;
interp1(x_nod,f(x_nod),3)
%%
M1 = 10; e1 = 1e-2; e2 = 1e-5;
M2 = ceil(M1*nthroot(e1/e2,4))
%%
f = @(x) sin(pi*x)+2*x^4-7*x^3;
integr_stima(-2,2,1e-4,f,'s')
%%
f = @(x) 5+x.^7; a = 0; b = 2;
gausslegendre_comp(a,b,f,1,2,'equi')
%%
syms b
h = 1/4;
e = -1/12*h^2*2*6*b
%%
f = @(t,y) -exp(y)+2*t; y0 = 3; tf = 10; tv = [0 tf]; h = 0.1;
[t,u] = Heun(f,10,y0,0.1);
u(2)
%%
clear
clc

mu = 1; eta = 5; sigma = @(x) 0.*x; f = @(x) 0.*x + 2; 
a = 0; b = 1; alpha = 0; beta = 3;
N = 9; h = 0.1;
[t,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
u(10)
%%
K1 = 1e4; h1 = 0.1; h2 = h1/10;
K2 = K1*(h1/h2)^2
%%
mu = 1; f = @(x,t) 3;
a = 0; b = 1; u_s = @(t) 0; u_d = @(t) 0; g_0 = @(x) 2*sin(pi*x);
T = 10; h = 1/2; delta_t = 1/8; theta = 1;
[u, x, t] = EqCalore_DiffFin_Theta(mu, f, a, b, u_s, u_d, g_0, T, h, delta_t,theta);
u(2,2)
%%
A = [-2 -8 0; 1 0 -1; 0 0 -1];
g = @(t) exp(-t/2).*[-3/2*(4*pi*sin(pi*t)+(4*pi^2-29)*cos(pi*t));2;1];
y0 = [-3 6 2]';
tf = 10; tv = [0 tf]; h = 0.1;
fun = @(t,y) A*y + g(t);
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(2,2),u(2,end)

l = max(eig(A));
fsolve(@(h) abs(1+h*l)-1,1)
[t,u] = Crank_Nicolson_sistemi_lineari(A,g,tv,y0,tf/h);
u(2,2),u(2,end)

y = @(t) [exp(-t/2).*(-3*cos(pi.*t)-6*pi*sin(pi.*t)); exp(-t/2).*(6*cos(pi.*t)); exp(-t/2).*2];
subplot(2,1,1)
plot(t,u,'LineWidth',2)
subplot(2,1,2)
plot(t,y(t),'LineWidth',2)
Eh = [];
h = [10 5 2.5 1.25].*1e-3;
for i = h
    [t,u] = Crank_Nicolson_sistemi_lineari(A,g,tv,y0,tf/i);
    Eh = [Eh norm(u(:,end)-y(tf))];
end
Eh,stimap_2(Eh,h);
figure()
loglog(h,Eh,h,h,'-',h,h.^2,':','LineWidth',2);
legend('Eh','H','H^2')
%%
clear
clc
f = @(x) abs(1-exp(abs(sin(x))));
n = 4; a = -2; b = 2;
x_dis = a:0.0001:b;
[y_dis, P, err] = poly_lagr(f,n,a,b,x_dis,'equi');
max(y_dis)
%%
clear
clc
f = @(x) x.*exp(sin(x));
a = 0; b = 10; x = a:0.0001:b;
M = 10; h = (b-a)/M;
x_nod = a:h:b;
[y_max,i] = max(abs(interp1(x_nod,f(x_nod),x)-f(x)));
y_max
x(i)
%%
f = @(x) exp(x);
a = -1; b = 1;
tol = 1e-2; 
integr_stima(a,b,tol,f,'pm')
%%
f = @(x) 5+x.^7;
a = 0; b = 3;
gausslegendre_comp(a,b,f,3,1,'equi')
%%
syms h;
f = @(x) x.^2.*(7.*x-3);
df = Jac(f);
ddf = Jac(df);
dddf = Jac(ddf);
e = -1/12*h^2*2*42
%%
f = @(t,y) -(1+t)*y;
y0 = 3;
h = 0.1;
[t,u] = Crank_Nicolson(f,10,y0,h);
u(2)
%%
clear
clc
mu = 1; eta = 0; sigma =@(x) 2+0.*x; f = @(x) 10*x;
a = 0; b = 1;
alpha = 0; beta = 1;
h = 0.1; N = 9;
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N);
u(6)
%%
h1 = 0.1; e1 = 4*1e-2; h2 = 0.05;
e2 = e1*(h2/h1)^2
%%
mu = 1; f = @(x,t) 0;
a = 0; b = 1; u_s = @(t) 0; u_d = @(t) 0; g_0 = @(x) 7*sin(pi*x); T = 10;
h = 0.5; delta_t = 0.2; theta = 1;
[u, x, t] = EqCalore_DiffFin_Theta(mu, f, a, b, u_s, u_d, g_0, T, h, delta_t,theta);
u(2,1/delta_t+1)
%%
a = -1; b = 1; f = @(x) 1./(1+9*x.^2);
[eq, phi, A] = poly_equation(10,a,b,'equi');
A
n = [6 8 10 12];
x_dis = a:0.0001:b; Eh = [];
for N = n
    [y_pol, P, err] = poly_lagr(f,N,a,b,x_dis,'CGL');
    Eh = [Eh max(abs(f(x_dis)-y_pol))];
end
Eh
%%
a = -1/4; b = -1/5; c = -1/8; d = 0;
y0 = [8 5]';
A = @(y) [a b*y(1);c*y(2) d];
g1 = @(t) 8*exp(-t/4)*((cos(pi*t))^2-pi*sin(pi*t));
g2 = @(t) 5*(exp(-t/4)*(cos(pi*t))^2-pi*sin(pi*t));
g = @(t) [g1(t) g2(t)]';
h = 0.05; tf = 10; tv = [0 tf];

fun = @(t,y) A(y)*y+g(t);
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(:,end)
[t,u] = BOH(A,g,tv,y0,tf/h);
u(:,end)