%% TDE
%% 2 ITINERE 2021/2022
% Test
%% 1
x = [0 0.5 1];
y = [2 1 1.5];
P = polyfit(x,y,2);
dP = polyder(P);
polyval(P,0.25)
polyval(dP,0.25)
%% 2
x = 0:0.25:1;
y = [3 0.5 1.5 -0.5 1];
x_dis = 0:0.001:1;
y_dis = interp1(x,y,x_dis);
plot(x_dis,y_dis)
min(y_dis)
%% 3
f = @(x) sin(pi.*x);
h = 2/4;
x_nod = 0:h:2;
interp1(x_nod,f(x_nod),1.6)
%% 4
x = 0 : 0.25 : 1;
y = [3 0.5 1.5 -0.5 1];
poly_scarto_IMQ(x,y,2)
%% 5
f = @(x) sqrt(2 + abs(x));
N = 3/1;
x_dis = -1:0.01:2;
trapcomp(-1,2,3,f)
%% 6
M1 = 20;
e1 = 1e-1;
e2 = 1e-3;
M2 = ceil(M1*nthroot(e1/e2,2))
%% 7
fun = @(t,y) -1/81*(1+sin(t))*y^2;
dfun = @(t,y) -2/81*(1+sin(t))*y;
y0 = 9;
h = 0.2;
[t,u] = eulero_indietro_newton(fun,dfun,5,y0,h);
plot(t,u)
u(5/h+1)
%% 8
mu = 1;
sigma = @(x) 2 + 0.*x;
f = @(x) 3 + x;
eta = 0;
[x,u] = miste(0,1,1,0,mu,eta,sigma,f,9,'1');
plot(x,u)
%% 9
mu = 1;
eta = -50;
sigma = @(x) 0.*x;
h = 1/10;
mu = mu*(1+abs(eta)*h/(2*mu));
f = @(x) 0.*x;
N = 9;
[x,u] = trasporto_diffusione(0,1,0,3,mu,eta,f,N);
plot(x,u)
u(2)
%% 10

%% Esercizio
%% 1
K = 1;
A = @(y) [-1 -(147/16+K*y(2)); 1 0];
g = @(t) [-9.*exp(-t./2).*(sin(3.*t).^2-1)-9/2.*exp(-t./4).*sin(3.*t);0];
tf = 5;
tv = [0 tf];
y0 = [-3/4; 3];
fun = @(t,y) A(y)*y + g(t);
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/1e-2);
sol1 = [u(:,2) u(:,end)]
[u2_min,i] = min(u(2,:));
sol2 = [u2_min t(i)]
y = @(t) 3*exp(-t/4).*[-1/4*cos(3*t)-3*sin(3*t);cos(3*t)];
h = [10 5 2.5 1.25].*1e-4;
Eh = [];
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh norm(u(:,end)-y(t(end)))];
end
Eh, stimap_2(Eh,h);

A_RK = [0 0; 2/3 0];
b = [1/4 3/4];
c = [0 2/3];
[t1,u1] = Runge_Kutta(A_RK,b,c,fun,tv,y0,tf/1e-2);
sol3 = [u1(:,2) u1(:,end)]
A = [-1 -(147/16); 1 0];
l = max(eig(A));
R_abs = @(h) abs(1+h*l+(h*l)^2/2)-1;
fsolve(R_abs,1)
fun = @(t,y) A*y;
tf = 430;
[t,u] = Heun_sistemi(fun,[0 tf],y0,tf/0.43);
plot(t,u(2,:),'LineWidth',2)
%% 1 APPELLO 2021/2022
% Test
%% 1
f = @(x) sin(x+sqrt(2));
[y,P] = poly_lagr(f,3,0,pi,2,'equi');
dP = polyder(P);
polyval(dP,2)
y
%% 2
f = @(x) exp(3.*x);
poly_stima(f,3,-1,1,'CGL')
%% 3
f = @(x) ((x./pi).^2).*(x>=0 & x<pi) + ((2-x./pi).^2).*(x>=pi & x< 2*pi);
x_nod = 2*pi/5.*(0:5);
x_dis = x_nod(1):0.001:x_nod(end);
y = cubicspline(x_nod,f(x_nod),x_nod(4))
%% 4
x = 0:4;
y = [2 2 0 1 0];
poly_minim_IMQ(1,x,y)
%% 5
f = @(x) exp(x);
a = -1;
b = 3;
I = trapcomp(-1,3,4,f)
%% 6
f = @(x) x.^3;
fun = @(x) 3*(f(3+x)+f(3-x))-(6^4)/4;
fsolve(fun,2)
%% 7
fun = @(t,y) -2*y;
y0 = 7;
[t,u] = teta_met(fun,[0 10],y0,10/0.1,3/4);
plot(t,u)
u(2)
%% Esercizio
z = @(t) exp(-t/2)*(2*cos(t)-7/2*sin(t))+40*exp(-t)*sin(t)^2;
f = @(y) [-2*y(1)-10*y(2)^2; y(1)];
g = @(t) [z(t);0];
fun = @(t,y) f(y) + g(t);
y0 = [2;0];
tf = 10;
tv = [0 tf];
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/1e-2);
sol1 = [u(2,2) u(2,end)];
x_es = @(t) 2.*exp(-t./2).*sin(t);
h = [10 5 2.5 1.25].*1e-4;
Eh = [];
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh max(abs(x_es(t)-u(2,:)))];
end
Eh
stimap_2(Eh,h);
J = @(t) [-2 -20*x_es(t); 1 0];
tf = 9;
t = 0:0.01:tf;
lambda = [];
for t_h = t
    lambda = [lambda max(eig(J(t_h)))];
end
h = 0:0.001:1;
R_abs = abs(1+h.*max(lambda));
figure()
plot(h,R_abs,'LineWidth',2)
yline(1)
R_abs = @(h) abs(1+h.*max(lambda))-1;
h_max = fsolve(R_abs,1)
tf = 10;
tv = [0 tf];
[t,u] = eulero_avanti_sistemi(fun,[0 tf],y0,tf/0.01);
figure()
plot(t,u(2,:),'LineWidth',2)
hold on
plot(t,x_es(t),'LineWidth',2)
legend('APPROX','ESATTO')

%% APPELLO 2 2021/2022
%% 1
f = @(x) sin(x+sqrt(5));
a = 0;
b = pi;
[y,P] = poly_lagr(f,3,a,b,1,'CGL');
y
polyval(polyder(P),1)
%% 2
f = @(x) sin(pi.*x);
poly_stima(f,3,-1,1,'equi')
%% 3
%% 4
x = 0:4;
y = [1 1 2 2 3];
polyval(polyfit(x,y,1),4)
%% 5
f = @(x) exp(x);
I = pmedcomp(-2,3,5,f)
%% 6
e1 = 0.08;
M2 = 50;
M1 = 25;
e2 = e1*(M1/M2)^2
%% 7
f = @(t,y) -3*y;
y0 = 4;
teta = 1;
h = 0.1;
[t,u] = teta_met(f,[0 10],y0,10/h,teta);
plot(t,u)
u(2)
%% 8
mu = 2;
f = @(x) 3*sin(x);
sigma = @(x) 1 + 0.*x;
h = 0.1;
a = 0; b = 1; alpha = 0; beta = sin(1);
[x_nodes,uh] = Dirichlet(a,b,alpha,beta,mu,sigma,f,9);
x_es = @(x) sin(x);
plot(x_nodes,uh,0:0.01:1,sin(0:0.01:1));
Eh = max(abs(x_es(x_nodes)-uh'))
%% Esercizio
A = [-1/6 -5; 1 0];
z = @(t) 0;
g = @(t) [z(t); 0];
fun = @(t,y) A*y + g(t);
h = 0.1;
tf = 50;
tv = [0 tf];
y0 = [-3/4 9]';
[t,u] = Heun_sistemi(fun,tv,y0,tf/0.01);
u(2,2)
u(2,end)
plot(t,abs(u(2,:)))
yline(1)
i = find(abs(u(2,:))>= 1,1,'last');
t_m = t(i+1)
x_t = @(t) 9.*exp(-t./12).*cos(sqrt(719)./12.*t);
h = [10 5 2.5 1.25].*1e-2;
err = [];
for i = h
    [t,u] = Heun_sistemi(fun,tv,y0,tf/i);
    err = [err max(abs(u(2,:)-x_t(t)))];
end
err
stimap_2(err,h)
h = 0:0.001:1;
l = max(eig(A));
R = @(h) abs(1 + h*l + (h*l)^2/2) - 1;
fsolve(R,1)
%% APPELLO 3 2021/2022
% Test
%% 1
f = @(x) sin(x+sqrt(7));
a = 0;
b = pi;
[~,P] = poly_lagr(f,3,a,b,pi,'equi');
dP = polyder(P);
ddP = polyder(dP);
polyval(dP,pi)
polyval(ddP,pi)
%% 2
f = @(x) exp(2.*(1-x));
n = 6;
poly_stima(f,n,0,2,'CGL');
x_dis = 0:0.001:2;
y_dis = poly_lagr(f,n,0,2,x_dis,'CGL');
k = 0:n;
a = 0; b = 2;
t = -cos(pi*k/n);
x_nod = ((b-a)/2)*t + (a+b)/2;
plot(x_nod,f(x_nod),'-o',x_dis,y_dis)
%% 3
f = @(x) x.^3;
[y,x,y0] = interp_tratti(f,-1,1,2,4,'equi',0.75);
y0
%% 4
n = 5;
x = 0:5;
y = [4 4 1 1 0 0];
poly_minim_IMQ(2,x,y)
%% 5
f = @(x) exp(x);
trapcomp(-3,1,4/1,f)
%% 6
M1 = 25; e1 = 0.08; e2 = 0.02;
M2 = M1*sqrt(e1/e2)
%% 7
f = @(t,y) -2*y;
y0 = 10;
teta = 3/4;
[t,u] = teta_met(f,[0 10],10,10/0.1,3/4);
plot(t,u);
u(2)
%% APPELLO 4 2020/2021
% Test
%% 1
x = 0:4;
y = [3 1 3 4 5];
polyval(polyfit(x,y,4),1.5)
%% 2 
f = @(x) exp(x);
n = 4;
poly_stima(f,4,0,2,'equi')
%% 3
f = @(x) x.^3 + sin(pi.*x);
y = poly_lagr(f,2,0,3,1.5,'equi')
%% 4
f = @(x) x.^2 + abs(x);
h = 2/4;
x_nod = -1:h:1;
x_dis = -1:0.001:1;
y = polyval(polyfit(x_nod,f(x_nod),2),x_dis);
plot(x_nod,f(x_nod),'o',x_dis,y,'LineWidth',2)

trapz(x_dis,y)
%% 5
M1 = 10;
M2 = 100;
e1 = 1e-1;
e2 = e1*(M1/M2)^2;
%% 6
syms t y;
f = -t^2*y;
y0 = 1;
EAI_sym(f,1,y0)
%% 7 
f = @(t,y) -(t+1)*y^3;
y0 = 1;
A = [0 0; 1/2 0];
b = [0 1];
c = [0 1/2];
[t,u] = Runge_Kutta(A,b,c,f,[0 10],y0,10/0.1);
plot(t,u)
u(2)
%% 8
A = [-2 1; 2 -3];
y0 = [1 3]';
tv = [0 1];
f = @(t,y) A*y;
[t,u] = Crank_Nicolson_sistemi(f,tv,y0,1/0.1);
plot(t,u)
u(:,end)
%% 9
A = [-2 -1; 1 0];
f = @(t,y) A*y;
y0 = [1 4]';
[t,u] = eulero_avanti_sistemi(f,[0 10],y0,10/0.1);
plot(t,u)
u(2,2)
%% 10
f = @(t,u,v) -u;
[t,u,v] = Leap_Frog(f,[0 10],0.1,5,0);
u(2)
u(3)
u(4)
%% APPELLO 1 2020
% Test
%% 11
x = linspace(0,1,50);
rng(1);
y = 2*x.^2 + 0.2*sin(100*pi*x) + 0.2*randn(1,50);
polyval(polyfit(x,y,2),0.5)
%% 12
f = @(x) (x+2).^2;
h = 10/5;
x_nod = 0:h:10;
interp1(x_nod,f(x_nod),1)
%% 13
f = @(x) sin(10.*x)-7.*x;
a = 0; b = 3;
df = Jac(f);
ddf = Jac(df);
toll = 1e-4;
max_df = max(abs(ddf(0:0.001:3)));
H = (sqrt(toll*8/max_df));
M = ceil(3/H)
%% 14
n = 6;
b = 1;
a = -1;
k = 0:n;
t = -cos(pi*k/n);
x_nod = ((b-a)/2)*t + (a+b)/2;
P = 1;
syms x
for i = 0: n
    P = P*(x-x_nod(i+1));
end
w = matlabFunction(simplify(P));
pmedcomp(-1,1,1,w)
%% 15
f = @(x) exp(x);
M = integr_stima(-2,2,1e-1,f,'t')
%% 17
f = @(t,y) -2*(t+1)*y^2;
tf = 10;
tv = [0 tf];
A = [0 0; 3/4 0];
b = [1/3 2/3];
c = [0 3/4];
[t,u] = Runge_Kutta(A,b,c,f,tv,3,tf/0.1);
u(2)
%% 18
f = @(t,y) -5*(t+1)*y;
tf = 1;
tv = [0 tf];
[t,u] = Crank_Nicolson(f,tf,8,0.2);
plot(t,u)
u(2)
%% 19
A = [-3 -4; 1 0];
y0 = [7 1]';
l = max(eig(A));
R_abs = @(h) abs(1+h*l);
fsolve(@(h)R_abs(h)-1,1)
%% 20
mu = 3;
eta = -50;
2*mu/abs(eta)
%% APPELLO 2 2020
%% Esercizio 
A = [-2 -6; 1 0];
g = @(t) [10;0];
y0 = [1 4]';
tf = 5;
tv = [0 tf];
fun = @(t,y) A*y + g(t);
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/0.1);
u(2,2)
u(2,end)
x = @(t) exp(-t)*(7/3*cos(sqrt(5)*t)+10/(3*sqrt(5))*sin(sqrt(5)*t))+5/3;
h = [10 5 2.5 1.25].*1e-3;
Eh = [];
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh abs(u(2,end)-x(t(end)))];
end
Eh
stimap_2(Eh,h);
[t,u] = Heun_sistemi(fun,tv,y0,tf/0.1);
u(2,2)
u(2,end)
l = max(eig(A));
R_abs = @(h) abs(1 + h*l + (h*l).^2./2);
h_max = fsolve(@(h)R_abs(h)-1,1)
%%
f = @(x) 3.*sin(x);
mu = 2;
eta = 0;
sigma = @(x) 1 + 0.*x;
[x_nodes,uh] = Dirichlet(0,1,0,sin(1),mu,eta,sigma,f,9,0.1);
u_ex = @(x) sin(x);
plot(x_nodes,uh,x_nodes,u_ex(x_nodes),'o')

max(abs(u_ex(x_nodes')-uh))
%%
t0=0;
tf=10;
h=1e-2;
Nh=round((tf-t0)/h);
y0=[2 0]';
[t_ea,u_ea]=eulero_avanti_sistemi(@es41,[t0 tf],y0,Nh);
x_ex=@(t) 2*exp(-t/2).*sin(t);
plot(t_ea,x_ex(t_ea),'b',t_ea,u_ea(1,:),'r',LineWidth=1.8)
grid on;
legend('Sol.esatta','Sol. EA',Location='best')
% valuto le approssimazioni u in t1 e tf
u_t1=u_ea(1,1)
u_tf=u_ea(1,tf)
% vedi il grafico (funzione es41)
% calcolo gli errori
passo=[1e-3 5e-4 2.5e-4 1.25e-4];
E=[];
e=[];
for i=1:length(passo)
    N=((tf-t0)/passo(i));
    [t,u]=eulero_avanti_sistemi(@es41,[t0 tf],y0,N);
    E=[E; max(abs(x_ex(t)-u(1,:)))];
  % e=[e; abs(x_h-u(1,:))];
end
E
figure(2);
loglog(passo,E,'r',passo,passo.^2,'k--',LineWidth=1.7)
grid on;
%% APPELLO 4 2020
% Esercizio
a = -1/4; b = -1/5; c = -1/8; d = 0;
y0 = [8 5]';
g1 = @(t) 8*exp(-t/4)*(cos(pi*t)^2-pi*sin(pi*t));
g2 = @(t) 5*(exp(-t/4)*cos(pi*t)^2-pi*sin(pi*t));
A = @(y) [a b*y(1); c*y(2) d];
g = @(t) [g1(t) g2(t)]';
fun = @(t,y) A(y)*y + g(t);

tf = 10;
tv = [0 tf];
h = 0.05;

[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
plot(t,u)
u(:,end)

[t,u] = appello2020(A,g, tv, y0, tf/h);
plot(t,u)
u(:,end)

u_ex = @(t) [8*exp(-t/4)*cos(pi*t); 5*cos(pi*t)];
h = [10 5 2.5 1.25].*1e-4;
Eh = [];
for i = h
    [t,u] = appello2020(A,g, tv, y0, tf/i);
    Eh = [Eh norm(u(:,end)-u_ex(tf))];
end
Eh
stimap_2(Eh,h)
%% 5 appello 2020
% Test
%% 13
a = 0; b = 5;
f = @(x) x + sin(pi*x+1);
n = 10;
x_nodes = a:(b-a)/n:b;
polyval(polyfit(x_nodes,f(x_nodes),1),6)
%% 14
f = @(x) exp(x);
a = 0; b = 2; e = 1e-2;
H = sqrt(e*8/f(b));
M = ceil((b-a)/H);
x_dis = a:0.001:b;
y_dis = interp1(a:(b-a)/M:b,f(a:(b-a)/M:b),x_dis);
Err = max(abs(f(x_dis)-y_dis))
%% 15
x = 1:5;
y = [1 1 0 1 2];
spline(x,y,4.5)
%% 16
syms b g;
fun = @(x) b*x.^3 + g*x.^2 + 1;
simplify(pmedcomp(-2,2,2,fun))
%% 17
a = -1; b = 1;
h = 1/2;
r = GaussLegendre_grado_r(a,b,h)
%% 18
f = @(x) 2.^x;
h = 1/4;
df = @(x) (-3*f(x)+4*f(x+h)-f(x+2*h))/(2*h);
df(0)
%% 19
f = @(t,y) 1+t*y;
y0 = 1;
tf = 5;
tv = [0 tv];
h = 0.1;
[t,u] = Crank_Nicolson(f,tf,y0,h);
u(2)
%% 20
mu = 1; sigma = @(x) (1 + 10*x); f = @(x) 5 + 0.*x;
h = 0.1;
N = 9;
[x_nodes,uh] = Dirichlet(0,1,0,0,mu,0,sigma,f,N,h);
plot(x_nodes,uh)
uh(6)
%% 22
f = @(t,x) 0.*x;
u_s = @(t) 0.*t;
u_d = u_s;
g_0 = @(x) 5*sin(pi*x);
[u,x,t] = EqCalore_DiffFin_Theta(0.5,f,0,1,u_s,u_d,g_0,1,0.5,0.1,0);
u(2,end)
%% Esercizio
m = 1; c = 2; k = @(t) 3*(2-exp(-t));
y0 = [-10 10]';
tf = 2;
tv = [0 tf];
z = @(t) 10*exp(-t)*cos(t)*(30*exp(-t)*(2-exp(-t))*cos(t)-2);
f = @(t,y) [-c/m*y(1)-k(t)/m*y(2)^2+z(t)/m; y(1)];

h = 0.05;
[t,u] = eulero_avanti_sistemi(f,tv,y0,tf/h);
u(:,end)

A = @(t,y) [-c/m -k(t)/m*y(2); 1 0];
g = @(t) [z(t)/m 0]';
f = @(t,y) A(y)*y + g(t);
[ t, u ] = appello5_2020( A, g, tv, y0, tf/h );
u(:,end)
y = @(t) [-10*exp(-t)*(cos(t)+sin(t)) 10*exp(-t)*cos(t)]';
h = [10 5 2.5 1.25].*1e-4;
Eh = [];
for i = h
    [ t, u ] = appello5_2020( A, g, tv, y0, tf/i );
    Eh = [Eh norm(u(:,end)-y(tf))];
end
format shortE
Eh
stimap_2(Eh,h)
%% KARIM
z = @(t) 0;
A = [-1/6 -5; 1 0];
g = @(t) z(t)/2;
fun = @(t,y) A*y + g(t);
y0 = [-3/4 9]';
tf = 50;
tv = [0 tf];
[t,u] = Heun_sistemi(fun,tv,y0,tf/1e-1);
pt2 = [u(2,2) u(2,end)]
figure()
plot(t,abs(u(2,:)),'-o','LineWidth',2)
yline(1)
g = find(abs(u(2,:))>=1,1,'last');
tm = t(g+1)

x = @(t) 9*exp(-t./12).*cos(sqrt(719)/12.*t);
h = [10 5 2.5 1.25].*1e-2;
Eh = [];
for i = h
    [t,u] = Heun_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh max(abs(u(2,:)-x(t)))];
end
p = log(Eh(end)/Eh(end-1))/log(h(end)/h(end-1))

l = max(eig(A));
R_abs = @(h) abs(1 + h.*l + (h.*l).^2/2);
phi = @(h) h - R_abs(h)+1;
h_max = ptofis(0.4,phi,1e5,1e-6);
h_max(end)
%h_max = fsolve(@(h) R_abs(h)-1, 1)
g = @(t) z(t)/2;
[t,u] = multipasso2(A,g,tv,y0,tf/1e-1);
pt6 = [u(2,2) u(2,3) u(2,end)]
%% 
f = @(y) [-2*y(1)-10*y(2)^2;y(1)];
z = @(t) exp(-t/2).*(2*cos(t)-7/2*sin(t))+40.*exp(-t).*sin(t).^2;
g = @(t) [z(t) 0]';
fun = @(t,y) f(y) + g(t);
y0 = [2 0]';
tf = 10;
tv = [0 tf];

[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/1e-2);
pt2 = [u(2,2) u(2,end)]
plot(t,abs(u(2,:)),'LineWidth',2)
yline(0.1,'LineWidth',2)
i = find(abs(u(2,:))>=0.1,1,'last');
t_m = t(i+1)

x = @(t) 2.*exp(-t./2).*sin(t);
h = 10.*2.^(0:-1:-3).*1e-4;
Eh = [];
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh max(abs(u(2,:)-x(t)))];
end
Eh
p = log(Eh(end)/Eh(end-1))/(log(h(end)/h(end-1)))
J = @(t) [-2 -20*x(t);1 0];
lambda = [];
for t = 0:0.001:tf
    lambda = [lambda max(eig(J(t)))];
end
lambda = max(lambda);
R_abs = @(h) abs(1 + h*lambda);
%% R_abs = 1 --> R_abs-1 = 0 --> phi(alpha) = alpha - f(alpha)
phi = @(h) h - R_abs(h) + 1;
h = 0:0.001:1;
plot(h,R_abs(h));
yline(1); ylim([0 2]);
succ = ptofis_plot(0.05,phi,1e5,1e-8,0,2);
h_max = succ(end)
[t,u] = multipasso(fun,tv,y0,tf/1e-2);
pt6 = [u(2,2) u(2,3) u(2,end)]

h = 10.*2.^(0:-1:-3).*1e-4;
Eh = [];
for i = h
    [t,u] = multipasso(fun,tv,y0,tf/i);
    Eh = [Eh max(abs(u(2,:)-x(t)))];
end
Eh
p = log(Eh(end)/Eh(end-1))/(log(h(end)/h(end-1)))
