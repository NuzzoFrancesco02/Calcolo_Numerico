%%
clear
clc
K = 1;
A = @(y) [-1 -(147/16+K*y(2)); 1 0];
g = @(t) [-9*exp(-t/2)*(sin(3*t).^2-1)-9/2*exp(-t/4)*sin(3*t),0]';
y0 = [-3/4 3]';
tf = 5;
tv = [0 tf];
A_RK = [0 0; 2/3 0];
b = [1/4 3/4];
c = [0 2/3];
h = 1e-2;
fun = @(t,y) A(y)*y+g(t);
[t,u] = Runge_Kutta(A_RK,b,c,fun,tv,y0,tf/h);
plot(t,u)
u(:,2)
u(:,end)
%%
clear
clc
A = diagonals([4 -2 -2],10);
g = @(t) cos(pi*t)*ones(10,1);
y0 = 7*ones(10,1);
fun = @(t,y) A*y + g(t);
tf = 10;
tv = [0 tf];
A_RK = [1/4 0; 1/2 1/4];
b = [1/2 1/2];
c = [1/4 3/4];
h = 0.1;
[t,u] = Runge_Kutta(A_RK,b,c,fun,tv,y0,tf/h);
plot(t,u,'LineWidth',2)

u(3,2)
u(3,end)
%%
clear
clc
A = [-4 -100; 1 0];
g = @(t) [200*cos(10*t) 0]';
fun = @(t,y) A*y + g(t);
y0 = [50 0]';
tf = 5;
tv = [0 tf];

A_RK = [0 0 0; 1/2 0 0; -1 2 0];
b = [1/6 2/3 1/6];
c = [0 1/2 1];
h = 1e-2;
[t,u] = Runge_Kutta(A_RK,b,c,fun,tv,y0,tf/h);
plot(t,u)
u(2,2)
u(2,end)
%%
A = diagonals([5/2 -2 -1/2],9);
g = @(t) (1 + 2*sin(pi*t))*ones(9,1);
y0 = 2*ones(9,1);
fun = @(t,y) A*y + g(t);
tf = 10;
tv = [0 tf];
A_RK = [1/4 0; 1/2 1/4];
b = [1/2 1/2];
c = [1/4 3/4];
h = 0.1;
[t,u] = Runge_Kutta(A_RK,b,c,fun,tv,y0,tf/h);
u(5,2)
u(5,end)
%%
f = @(t,y) -(t+1)*y^2;
y0 = 9;
h = 0.1;
tf = 10;
tv = [0 tf];
A = [0 0; 1 0];
b = [1/2 1/2];
c = [0 1];
[t,u] = Runge_Kutta(A,b,c,f,tv,y0,tf/h);
u(2)
%%
A = diagonals([3 -2 -1],9);
g = @(t) exp(sin(pi*t))*ones(9,1);
y0 = 4*ones(9,1);
fun = @(t,y) A*y + g(t);
A_RK = [1/3 0; 1/3 1/3];
b = [1/2 1/2];
c = [1/3 2/3];
tf = 10;
tv = [0 tf];
[t,u] = Runge_Kutta(A_RK,b,c,fun,tv,y0,tf/0.1);
plot(t,u)
u(5,2), u(5,end)
%%
f = @(t,y) -(t+1)*y^3;
h = 0.1;
tf = 10;
tv = [0 tf];
y0 = 1;
A = [0 0; 1/2 0];
b = [0 1];
c = [0 1/2];
[t,u] = Runge_Kutta(A,b,c,f,tv,y0,tf/h);
u(2)
%%
A = [0 0; 3/4 0];
b = [1/3 2/3];
c = [0 3/4];
tf = 10;
tv = [0 tf];
y0 = 3;
fun = @(t,y) -2*(t+1)*y^2;
h = 0.1;
[t,u] = Runge_Kutta(A,b,c,fun,tv,y0,tf/h);
format long
u(2)
%%
f = @(t,y) -5*(t+1)*y;
y0 = 8;
theta = 1/2;
[t,u] = Crank_Nicolson(f,10,y0,0.2);
u(2)
format short
%%
A = [-2 -6;1 0];
f = @(t) 10;
g = @(t) [f(t) 0]';
y0 = [1 4]';
tf = 5; tv = [0 tf];
fun = @(t,y) A*y + g(t);
h = 0.1;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(2,2), u(2,end)
x = @(t) exp(-t)*(7/3*cos(sqrt(5)*t)+10/(3*sqrt(5))*sin(sqrt(5)*t))+5/3;
h = [10 5 2.5 1.25].*1e-3;
Eh = [];
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh abs(u(end)-x(t(end)))];
end
Eh
stimap_2(Eh,h);
[t,u] = Heun_sistemi(fun,tv,y0,tf/0.1);
u(2,2)
u(2,end)
l = max(eig(A));
h = 0:0.01:1;
R = abs(1 + h.*l + (h.*l).^2/2);
plot(h,R,'LineWidth',2)
yline(1,'LineWidth',2)
fsolve(@(h) abs(1 + h.*l + (h.*l).^2/2)-1,1)
%% 
mu = 3;
h = 0.2;
delta_t = h^2/(2*mu)
%%
options = optimset('Display','off');
A = [-1 -(147/16);1 0];
tf = 5;
l = max(eig(A));
R_abs = @(h) abs(1 + h*l + (h*l)^2/2)-1;
fsolve(R_abs,1,options)
h = 0:0.001:1;
R_abs = abs(1 + h.*l + (h.*l).^2./2);
plot(h,R_abs,'linewidth',2);
yline(1,'LineWidth',2)
%%
A = [-1 -(147/16);1 0];
tf = 5;
fun = @(t,y) A*y;
teta = 1/2;
h = 1e-2;

[t,u1] = teta_met(fun,[0 tf],[-3/4; 3],tf/h,teta);
plot(t,u1)
u1(:,2)
u1(:,end)
[t,u2] = Crank_Nicolson_sistemi(fun,[0 tf],[-3/4; 3],tf/h);
u2(:,2)
u2(:,end)
%%
mu = 1;
eta = 0;
sigma = @(x) 0.*x + pi^2;
f = @(x) pi^2.*(x+2*sin(pi.*x));
h = 0.1;

a = 0; b = 1; alpha = 1 + pi; beta = 1 - pi; N = (b-a)/h-1;
[x_nodes, u] = miste(a,b,alpha,beta,mu,eta,sigma,f,N,'2','N-N');
u(end)
%%
mu = 1; eta = 0; sigma = @(x) 1 + 19.*x; f = @(x) 0.*x + 5; a = 0; b = 1; alpha = 0; beta = 0;
[x,u] = miste(a,b,alpha,beta,mu,eta,sigma,f,N,'1','N-N');
u(6)
%%
mu = 1; a = 0; b = 1; alpha = 1; gamma = 0;
eta = 0;
sigma = @(x) 0.*x + 2;
f = @(x) 3 + x; 
h = 0.1;
N = 9;
[x, u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,'2','D-N');
u(end)
%%
mu = 1;
eta = -100;
sigma = @(x) 0.*x;
f = @(x) 0.*x;
a = -1; b = 1; alpha = 7; beta = 0;
h = 0.2;
Peh = abs(eta)*h/(2*mu);
mu = mu*(1+Peh);
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,(b-a)/h-1,h);
u(2)
%%
fdt = tf([1],conv([1 1],[1 2 1]));
bode(fdt)
%%
f = @(t,y) -y+2*t;
y0 = 1;
EA_sym(f,3,y0)
%%

x = [0 0.5 1]; y = [2 1 1.5];
P = polyfit(x,y,2);
dP = polyder(P);
polyval(P,0.25)
polyval(dP,0.25)
%%
x = 0:0.25:1; y = [3 0.5 1.5 -0.5 1];
min(y)
%%
f = @(x) sin(pi*x); n = 4; h = 2/n; x_nod = 0:h:2;
interp1(x_nod,f(x_nod),1.6)
%%
n = 4; x = 0:0.25:1; y = [3 0.5 1.5 -0.5 1];
poly_scarto_IMQ(x,y,2)
%%
f = @(x) sqrt(2+abs(x));
a = -1; b = 2;
N = 3;
trapcomp(a,b,N,f)
%%
M1 = 20; e1 = 1e-1; e2 = 1e-3;
M2 = ceil(M1*nthroot(e1/e2,2))
%%
f = @(t,y) -(1+sin(t))*y^2/81;
y0 = 9;
df = @(t,y) -2*(1+sin(t))*y/81;
tf = 5; tv = [0 tf]; h = 0.2;
[t,u] = eulero_indietro_newton(f,df,tf,y0,h);
u(end)
%%
mu = 1; eta = 0; sigma = @(x) 0.*x+2; f = @(x) 3+x;
h = 0.1; a = 0; b = 1; alpha = 1; gamma = 0; N = 9;
[x,u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,N,'1','D-N','',h);
u(end)
%%
mu = 1; eta = -50; sigma = @(x) 0.*x; f = @(x) 0.*x; a = 0; b = 1; alpha = 0; beta = 3;
h = 0.1; N = 9;
Peh = abs(eta)*h/(2*mu);
mu = mu*(1+Peh);
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
u(2)
%%
h = 0.2; mu = 3;
h^2/(2*mu)
%%
K = 1;
A = @(y) [-1 -(147/16+K*y(2)); 1 0];
g = @(t) [-9*exp(-t/2)*(sin(3*t)^2-1)-9/2*exp(-t/4)*sin(3*t) 0]';
fun = @(t,y) A(y)*y + g(t);
tf = 5; tv = [0 tf];  y0 = [-3/4 3]';
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(:,2),u(:,end)
[u_min, i] = min(u(2,:));
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
A_RK = [0 0; 2/3 0]; b = [1/4 3/4]; c = [0 2/3];
h = 1e-2;
[t,u] = Runge_Kutta(A_RK,b,c,fun,tv,y0,tf/h);
u(:,2),u(:,end)
A = [-1 -147/16; 1 0];
l = max(eig(A));
fsolve(@(h) abs(1 + h*l + (h*l)^2/2)-1,1)
fun = @(t,y) A*y;
[t,u] = Crank_Nicolson_sistemi(fun,tv,y0,tf/1e-2);
u(:,2),u(:,end)

%%
f = @(x) sin(x+sqrt(2));
[y,P] = poly_lagr(f,3,0,pi,2,'equi'); y
polyval(polyder(P),2)
%%
f = @(x) exp(3*x);
poly_stima(f,3,-1,1,'CGL')
%%
f = @(x) ((x/pi).^2).*(x>=0 & x<pi) + ((2-x/pi).^2).*(x>=pi & x<2*pi);
n = 5; i = 0:n;
x = 2*pi/n.*i;
cubicspline(x,f(x),x(4))
%%
x = 0:4; y = [2 2 0 1 0];
poly_minim_IMQ(1,x,y)
%%
f = @(x) exp(x); a = -1; b = 3; N = 4; trapcomp(a,b,N,f)
%%
f = @(x) x.^3; 
y0 = 6^4/4;
fsolve(@(c) 3*(f(3+c)+f(3-c))-y0,1)
%%
fun = @(t,y) -2*y; y0 = 7; teta = 3/4; tf = 10; tv = [0 tf]; h = 0.1;
[t,u] = teta_met(fun,tv,y0,tf/h,teta);
u(2)
%%
h1 = 1e-2; e1= 4*1e-3; e2 = 1e-3;
h2 = h1*nthroot(e2/e1,2)
%%
mu = 1; eta = 0; sigma = @(x) (1+19*x); f = @(x) 0.*x + 5;
a = 0; b = 1; alpha = 0; beta = 0;
h = 0.1; N = 9;
[x,u] = miste(a,b,alpha,beta,mu,eta,sigma,f,N,'1','N-N','',h);
u(6)
%%
mu = 1; eta = 0; sigma = @(x) 0.*x; f = @(x) 0.*x+6;
a = 0; b = 1; b_b = 3; alpha = 0; beta = 0; h = 0.1; N = 9;
[x,u] = Robin(a,b,alpha,gamma,b_b,mu,eta,sigma,f,N,h);
u(end)
%%
clear
clc

z = @(t) exp(-t/2)*(2*cos(t)-7/2*sin(t))+40*exp(-t)*sin(t)^2;
f = @(y) [-2*y(1)-10*y(2)^2 y(1)]';
g = @(t) [z(t) 0]';
fun = @(t,y) f(y)+g(t);
tf = 10; tv = [0 tf];
y0 = [2 0]'; h = 1e-2;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
u(2,2),u(2,end)
plot(t,abs(u(2,:)),'o'),yline(0.1)
i = find(abs(u(2,:))>=0.1,1,'last');
t(i+1)
%%
clear
clc
z = @(t) exp(-t/2)*(2*cos(t)-7/2*sin(t))+40*exp(-t)*sin(t)^2;
f = @(y) [-2*y(1)-10*y(2)^2 y(1)]';
g = @(t) [z(t) 0]';
fun = @(t,y) f(y)+g(t);
tf = 10; tv = [0 tf];
y0 = [2 0]';

h = [10 5 2.5 1.25].*1e-4; Eh = [];
x = @(t) 2.*exp(-t./2).*sin(t);
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);  
    Eh = [Eh max(abs(u(2,:)-x(t)))];
end
Eh, stimap_2(Eh,h);

[t,u] = multipasso(fun,tv,y0,tf/1e-2);
u(2,2),u(2,3),u(2,end)
h = [10 5 2.5 1.25].*1e-4; Eh = [];
x = @(t) 2.*exp(-t./2).*sin(t);
for i = h
    [t,u] = multipasso(fun,tv,y0,tf/i);  
    Eh = [Eh max(abs(u(2,:)-x(t)))];
end
Eh, stimap_2(Eh,h);
%%
x = @(t) 2.*exp(-t./2).*sin(t);
J = @(t) [-2 1;-20*x(t) 0];
t = 0:0.001:tf; l =[];
for i = 1:length(t)
    l = [l max(eig(J(t(i))))];
end
l = max(l);
fsolve(@(h) abs(1+h*l)-1,0.1)
h = 0:0.001:1;
plot(h,abs(1+h.*l)),yline(1)
%%
f = @(x) sin(x+sqrt(5));
[y,P] = poly_lagr(f,3,0,pi,1,'CGL'); y
polyval(polyder(P),1)
%%
f = @(x) sin(pi*x); a = -1; b = 1; poly_stima(f,3,a,b,'equi')
%%
0
%%
x = 0:4; y = [1 1 2 2 3]; polyval(polyfit(x,y,1),4)
%%
f = @(x) exp(x); a = -2; b = 3; N = 5; pmedcomp(a,b,N,f)
%%
M1 = 25; e1= 0.08; M2 = 50; e2 = e1*(M1/M2)^2
%%
f = @(t,y) -3*y; teta = 1; y0 = 4; h = 0.1; tf = 10; tv = [0 tf];
[t,u] = teta_met(f,tv,y0,tf/h,teta); u(2)
[t,u] = eulero_indietro(f,tf,y0,h); u(2)
%%
mu = 2; eta = 0; sigma = @(x) 0.*x+1; f = @(x) 3*sin(x);
a = 0; b = 1; alpha = 0; beta = sin(1); h = 0.1; N = 9;
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
u_ex = @(x) sin(x);
max(abs(u'-u_ex(x)))
%%
mu = 1; eta = 100; sigma = @(x) 0.*x; f = sigma;
a = 0; b = 1; alpha = 5; beta = 0; h = 0.1; N = 9;
Peh = abs(eta)*h/(2*mu);
mu= mu*(1+Peh)
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
u(10)
%%
mu = 1; eta = 0; sigma = @(x) 0.*x; f = @(x) 0.*x+4; a = 0; b = 1; b_b = 3; alpha = 0; beta = 2;
N = 9; h = 0.1;
[x,u] = Robin(a,b,alpha,beta,b_b,mu,eta,sigma,f,N,h);
u(end)
%%
A = [-1/6 -5; 1 0]; 
z = @(t) 0;
g = @(t) [z(t)/2 0]';
y0 = [-3/4 9]';
fun = @(t,y) A*y+g(t);
tf = 50; tv = [0 tf]; h = 1e-1; 
[t,u] = Heun_sistemi(fun,tv,y0,tf/h);
u(2,2),u(2,end)
plot(t,abs(u(2,:)),'-o'),yline(1)
i = find(abs(u(2,:))>=1,1,'last');
t(i+1)

x = @(t) 9.*exp(-t./12).*cos(sqrt(719)./12.*t);
h = [10 5 2.5 1.25].*1e-2; Eh = [];
for i = h
    [t,u] = Heun_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh max(abs(u(2,:)-x(t)))];
end
Eh, stimap_2(Eh,h);
l = max(eig(A));
h = 0:0.001:1;
plot(h,abs(1+h.*l+(h.*l).^2/2)),yline(1)
fsolve(@(h) abs(1+h*l+(h*l)^2/2)-1,0.3)
%%
h = 1e-1;
[t,u] = multipasso2(A,g,tv,y0,tf/h);
u(2,2),u(2,3),u(2,end)
h = [10 5 2.5 1.25].*1e-2; Eh = [];
for i = h
    [t,u] = multipasso2(A,g,tv,y0,tf/i);
    Eh = [Eh max(abs(u(2,:)-x(t)))];
end
Eh, stimap_2(Eh,h);
%%
syms h y0 l;
s = 4;
esp = 0:s;
R = 1 + h*l + (h*l)^2/2 + (h*l)^3/6 + (h*l)^4/factorial(4);
max = solve(abs(R)-1,h,"Real",true);
max = double(max*l);
abs(max)