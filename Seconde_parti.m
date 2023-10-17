%% 1
x = [0 0.5 1];
y = [2 1 1.5];
x_dis = linspace(0,1,1000);
y_dis = polyval(polyfit(x,y,2),0.25)
%%
x = [ 0 0.25 0.5 0.75 1];
y = [3 0.5 1.5 -0.5 1];
x_dis = 0:0.01:1;
min(interp1(x,y,x_dis))
%%
f = @(x) sin(pi.*x);
h = 2/4;
x = 0:h:2;
interp1(x,f(x),1.6)
%%
n = 4;
x = 0:0.25:1;
y = [3 0.5 1.5 -0.5 1];
p = polyval(polyfit(x,y,2),x);
sum((p-y).^2)
%%
f = @(x) sin(x+sqrt(2));
a = 0;
b = pi;
lagr_polin(3,f,a,b,2,'equi')
%%
f = @(x) exp(3.*x);
CGL_stim(f,-1,1,3)
%%
f = @(x) ((x./pi).^2).*(x<pi & x >= 0) + ((2-x./pi).^2).*(x>=pi & x < 2*pi);
i = 0:5;
x = 2*pi.*i./5;
cubicspline(x,f(x),x(4))
%%
x = [0 1 2 3 4];
y = [2 2 0 1 0];
P = polyequation(x,y)
%%
x = [0 1 2 3 4];
y = [2 2 0 1 0];
minim_IMQ(1,x,y)
%%
f = @(x) sin(x+sqrt(5));
a = 0;
b = pi;
y = lagr_polin(3,f,a,b,1,'CGL')
%%
x = 0:4;
y = [1 1 2 2 3];
polyval(polyfit(x,y,1),4)
%%
f = @(x) sin(x+sqrt(7));
lagr_polin(3,f,0,pi,pi,'equi')
%%
f = @(x) x.^3;
a = -1;
b = 1;
h = (b-a)/4;
x_nod = [0 0.5 1];
y = polyval(polyfit(x_nod,f(x_nod),2),0.75)
%%
x = [0 1 2 3 4 5]';
y = [4 4 1 1 0 0]';
minim_IMQ(2,x,y)
%%
f = @(x) exp(3.*x);
x = -1:0.001:1;
[y,err] = lagr_polin(3,f,-1,1,x,'CGL');
max(err)
plot(x,f(x),x,y)
%%
x = [0 1 2 3 4];
y = [2 2 0 1 0];
n = 4;
minim_IMQ(1,x,y)
%%
f = @(x) exp(2.*(1-x))
poly_stima(f,6,0,2,'lagr')
%%
clear
clc
a = -1;
b = 1;
f = @(x) x.^3;
x = a:0.001:b;
[y,y0,x] = interp_quadr(f,a,b,4,'equi',0.75);
plot(x,f(x),x,y)
y0
%%
x = 0:5;
y = [4 4 1 1 0 0];
poly_minim_IMQ(2,x,y)
%%
f = @(x) cos(pi.*x);
sti = poly_stima(f,4,-1,1,'equi')
%%
x = 0:4;
y = [5 5 0 0 0];
poly_minim_IMQ(1,x,y)
%%
f = @(x) exp(2.*x + sin(pi.*x));
x = -1:0.001:1;
y = poly_lagr(f,5,-1,1,x,'equi');
[err_max,i] = max(abs(f(x)-y));
err_max, x(i)
%%
x = 0:0.25:1;
y = [1 0.5 1.5 -0.25 1];
pt = 0:0.0001:1;
fx = interp1(x,y,pt);
max(fx)
%%
f = @(x) 2.*abs(sin(pi.*x));
[y,x,y0] = interp_tratti(f,0,4,2,4,'equi',1.75);
plot(x,f(x),x,y)
%%
f = @(x) exp(2.*x)
poly_stima(f,5,-1,1,'equi')
%%
giorni = 1:8;
casi = [679 776 882 794 932 808 480 907];
plot(giorni,casi,'o-','LineWidth',2,'MarkerSize',10)
hold on;
plot(1:0.01:13,polyval(polyfit(giorni,casi,2),1:0.01:13),'LineWidth',2)
polyval(polyfit(giorni,casi,2),10)
%%
f = @(x) exp(x./2);
[y,x,y0] = interp_quadr(f,0,3,3,'equi',0.75);
plot(x,f(x),x,y,'LineWidth',2)
y0
%%
f = @(x) sqrt(x + sin(2*pi.*x) + 3);
poly_lagr(f,5,-1,1,0.4,'equi')
%%
f = @(x) exp(x./6);
[y,x,y0]=interp_tratti(f,0,3,4,3,'equi',1.5);
y0
plot(x,f(x),x,y,'LineWidth',2)
%%
f = @(x) exp(x).*(1+sin(pi.*x));
x = -1:0.0001:1;
y = poly_lagr(f,5,-1,1,x,'equi');
max(abs(f(x)-y))
%%
T = [6.0 5.5 6.0 4.5 4.0 3.5 4.0 4.5 5.0];
giorni = 1:9;
Y = polyval(polyfit(giorni,T,2),1:0.001:13);
plot(giorni,T,'LineWidth',2,'MarkerSize',10);
hold on;
plot(1:0.001:13,Y,'LineWidth',2)
polyval(polyfit(giorni,T,2),10)
%%
x = 0:0.25:1;
y = [3 0.5 1.5 -0.5 1];
p2 = polyval(polyfit(x,y,2),x);
scarto = sum((p2-y).^2)
%%
f = @(x) sqrt(2+abs(x));
N = 3*1;
trapcomp(-1,2,3,f)
%%
f = @(x) exp(x);
a = -1;
b = 3;
x = a:0.001:b;
x_nod = a:1:b;
plot(x,f(x),x,interp1(x_nod,f(x_nod),x))
fun = @(x) interp1(x_nod,f(x_nod),x)
trapcomp(a,b,4,fun)
%%
e2 = ((25^2)/(50^2))*0.08
%%
f = @(x) abs(sin(pi.*x));
H = 1;
N = 100;
simpcomp(0,100,100,f)
(100*0.01)/(20^2)
%%
eq = poly_equation([0 0.5 2]);


%%
log(1e-3)/log(10/100)
%%
f = @(x,y) 2.^(x+3.*y);
trap_doppio(0,1,0,1,f)
%%
n = 149;
x = linspace(0,1,150);
rng(1);
y = 3 * x.^2 + 0.3 * sin(100*pi*x)+0.3*randn(1,150);
polyfit(x,y,2)
p = polyval(polyfit(x,y,2),1.5)
%%
exp(x)
a = -1;
b = 2;
N = 1;
h = (b-a)/N;
x_nod = a:h:b;

f = @(x) polyval(polyfit(x_nod,exp(x_nod),N),x);
trapcomp(a,b,N,f)
%% 
e2 = (1e-1)*(10/100)^2
%%
f = @(x,y) exp((2.*x) + y);
integr_doppio(2,2,f,'pm_ret',0,1,0,1)
%%
a = 0;
b = 1;
n = 3;
[~,phi] = poly_equation(a:(b-a)/n:b)
pretty(expand(phi))
%%
f = @(x) 1 + x.^5 + abs(x);
a = -1;
b = 1;
n = 7;
x_nod = a:(b-a)/n:b;
s = @(x) cubicspline(x_nod,f(x_nod),x);
trapcomp(a,b,1,s)
%%
e2 = 0.1*(20/40)^2
%%
f = @(x,y) exp(2.*x+ y);
integr_doppio([0 1 0],[0 0 1],f,'pm_tri')
%%
a = -1;
b = 1;
f = @(x) exp(x).*(1+sin(pi.*x));
n = 5;
x = a:0.001:b;
y = poly_lagr(f,n,a,b,x,'equi');
max(abs(f(x)-y))
%%
T = [6 5.5 6 4.5 4 3.5 4 4.5 5];
giorni = 1:length(T);
plot(giorni,T,'o','MarkerSize',10,'LineWidth',2);
x = 1:0.001:15;
y = polyval(polyfit(giorni,T,2),x);
hold on
grid on
plot(x,y,'LineWidth',2);

polyval(polyfit(giorni,T,2),10)
%%
f = @(x) sqrt(1+x);
M = 4;
h = (4)/4;
x = -1:h:3;
trapcomp(-1,3,4,f)
%%
f = @(x,y) exp(x+3.*y);
integr_doppio([],[],f,'GL_1',0,1,0,1)
%%
f = @(x) x.^2 + abs(x);
a = -1; 
b = 1;
h = (b-a)/4;
x_nod = a:h:b;
x = a:0.001:b;
fp = @(x) polyval(polyfit(x_nod,f(x_nod),2),x);
plot(x,f(x),x,fp(x))
trapcomp(a,b,1000,fp)
%%
a = 0;
b = 10;
f = @(x) (x+2).^2;
h = (b-a)/5;
x_nod = a:h:b;
interp1(x_nod,f(x_nod),1)
%%
f = @(x) sin(10.*x) - 7.*x;
a = 0;
b = 3;
%%
n = 6;
i = 0:n;
x_nod = -cos(pi.*i./n);
        syms y;
        omega = 1;
        for i = 1 : n + 1
            omega = omega*(y-x_nod(i));
        end
        omega = matlabFunction(simplify(expand(omega)));
pmedcomp(-1,1,6,omega)
%%
f = @(x) exp(x);
a = -2;
b = 2;
%%
f = @(x) cos(3.*x+sqrt(5));
poly_stima(f,4,0,1,'equi')
%%
f = @(x) exp(x);
x_nod = -2:1:2;
[~,~,y0]=interp_tratti(f,-2,2,2,4,'equi',1.5)
%%
x = 0:10;
f = @(x) exp(x./10)+0.1*sin(pi.*x+sqrt(2));
polyval(polyfit(x,f(x),2),11)
%%
a = -2;
b = 2;
f = @(x) abs(1-exp(abs(sin(x))));
x = a:0.001:b;
y = poly_lagr(f,4,a,b,x,'equi');
plot(x,f(x),x,y)
max(y)
%%
f = @(x) x.*exp(sin(x));
x_nod = linspace(0,10,11);
x = 0:0.001:10;
y = interp1(x_nod,f(x_nod),x);
[e,i] = max(abs(f(x)-y))
x(i)
plot(x,f(x),x,y)
%%
f = @(x) exp(3.*x);
poly_stima(f,3,-1,1,'CGL')
%%
f = @(x) x.^3
sol = @(c) 3.*(f(3+c))+3.*f(3-c)-6^4/4;
fsolve(sol,1.5)
%%
f = @(x) sin(x);
fp = matlabFunction(poly_trigo(f,3));
%[y,x] = interp_tratti(f,-1,1,2,10000,'CGL',10);
x = -pi:0.001:pi;
subplot(2,1,1)
plot(x,f(x))
subplot(2,1,2)
plot(x,fp(x))
%%
f = @(x) exp(2*x+sin(pi.*x));
x = -1:0.0001:1;
y = poly_lagr(f,5,-1,1,x,'equi');
[e_max, i] = max(abs(f(x)-y));
e_max
x(i)
%%
x = 0:0.25:1;
y = [1 0.5 1.5 -0.25 1];
in = 0:0.001:1;
yp = interp1(x,y,in);
plot(in,yp)
hold on
plot(x,y,'o')
%%
f = @(x) 2.*abs(sin(pi.*x));
[y,x,x0] = interp_tratti(f,0,4,2,4,'equi',1.75);
plot(x,f(x),x,y)
x0
%%
x = [0 0.5 2];
[eq,phi]= poly_equation(x);
pretty(phi)
phi0 = matlabFunction(phi(1));
[y,x]=interp_tratti(phi0,0,2,2,1,'equi',10);
plot(x,phi0(x),x,y)
simpcomp(0,2,1,phi0)
%%
f = @(x,y) 2.^(x+3.*y);
integr_doppio([],[],f,'t',0,1,0,1)
%%
f = @(x) sin(x+sqrt(2));
a = 0;
b = pi;
df = Jac(f);
y = poly_lagr(f,3,a,b,2,'equi')
y = poly_lagr(df,3,a,b,2,'equi')
%%
f = @(x) exp(3.*x);
poly_stima(f,3,-1,1,'CGL');
%%
x = 0:4;
y = [2 2 0 1 0];
poly_minim_IMQ(1,x,y)
%%
f = @(x) exp(x);
x = -1:1:3;

y = interp1(x,f(x),-1:0.001:3);
trapcomp(-1,3,4,f)
%%
f = @(x) sin(x+sqrt(5));
poly_lagr(f,3,0,pi,1,'CGL')
df = Jac(f);
poly_lagr(df,3,0,pi,1,'CGL')
%%
f = @(x) sin(pi.*x);
poly_stima(f,3,-1,1,'equi')
%%
f = @(x) sin(x+sqrt(7));
x = pi;
[~,P] = poly_lagr(f,3,0,pi,x,'equi');
dP = polyder(P);
ddP = polyder(dP);
polyval(dP,pi)
polyval(ddP,pi)
%%
f = @(x) exp(2.*(1-x));
poly_stima(f,6,0,2,'equi')
%%
h = 2/4;
x = -1:h:1;
f = @(x) x.^3;
[y,x,y0]=interp_tratti(f,-1,1,2,4,'equi',0.75);
y0
%%
x = 0:5;
y = [4 4 1 1 0 0];
poly_minim_IMQ(2,x,y)
%%
f = @(x) cos(x);
[y,P] = poly_lagr(f,3,0,pi,pi,'equi');
polyval(polyder(P),pi)
polyval(polyder(polyder(P)),pi)
%%
f = @(x) cos(pi.*x);
poly_stima(f,4,-1,1,'equi')
%%
x = [0:4];
y = [5 5 0 0 0];
poly_minim_IMQ(1,x,y)
%%
f = @(x) abs(sin(pi.*x));
simpcomp(0,100,100,f)
%%
f = @(t,y) -sqrt((10.*y)./(10+t));
y0 = 9;
h = 0.2;
[t,u] = Heun(f,4,y0,h);
u(end)
%%
fun = @(t,y) -2*y+0*t;
[t,u] = teta_met(fun,[0 0.1],7,1,3/4);
u
%% 
Af = diagonals([4 -2 -2],10);

f = @(t,y) Af*y + cos(pi*t);
tf = 10;

A = [1/4 0; 1/2 1/4];
b = [1/2 1/2];
c = [1/4 3/4]';
%%
[t,u] = Runge_Kutta_2(f,tf,7*ones(1,10),0.1);
u(1,3)
u(7,3)

%%
A = [0];
b = [1];
c = [0]';
[~] = Runge_Kutta_gen(A,b,c);
