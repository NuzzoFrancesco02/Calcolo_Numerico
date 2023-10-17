%%
z = @(t) 0;
A = @(y) [-1/12 -4; 1 0];
g = @(t) [z(t)/3 0]';
fun = @(t,y) A(y)*y+g(t);
[x,u] = eulero_avanti_sistemi(@prova,[0 100],[-5/12 10]',100/(1e-2));
f = @(t) 10.*exp(-t./24).*cos(sqrt(2303/576).*t);
plot(x,abs(u(2,:)))
yline(2)
%%
close all
h = [1e-3 5*1e-4 2.5*1e-4 1.25*1e-4];
err = [];
for i = h
    [x,u] = eulero_avanti_sistemi(fun,[0 100],[-5/12 10]',100/i);
    err = [err max(abs(f(x)-u(2,:)))];
end
p = stimap_2(err,h);
h = 0:0.001:1;
l = max(eig([-1/12 -4; 1 0]));
z = 1 + h.*l + ((h.*l).^2)./2;
R = (real(z).^2+imag(z).^2);
plot(h,R)
yline(1)
R_abs = @(h) abs((real(1 + h.*l + ((h.*l).^2)./2).^2+imag(1 + h.*l + ((h.*l).^2)./2).^2))-1;
fsolve(R_abs,0.3)
[x,u] = esame([-1/12 -4; 1 0],g,[0 100],[-5/12 10]',100/1e-1);
%plot(t,u(2,:))
u(2,2)
u(2,3)
u(2,end)
h = [1e-3 5*1e-4 2.5*1e-4 1.25*1e-4];
err = [];
for i = h
    [x,u] = esame([-1/12 -4; 1 0],g,[0 100],[-5/12 10]',100/i);
    err = [err max(abs(f(x)-u(2,:)))];
end
p = stimap_2(err,h);
p(end)
%%
f = @(x) exp(x).*(1+sin(pi.*x));
n = 5;
x = -1:0.001:1;
y = poly_lagr(f,n,-1,1,x,'equi');
err_max = max(abs(y-f(x)))
%%
T = [6 5.5 6 4.5 4 3.5 4 4.5 5];
g = 1:9;
T_10 = polyval(polyfit(g,T,2),10)
%%
syms g x;
der = (g-1)*x^(g-1);
n = 1;
b = 1;
a = 0;
(abs(der))*(((b-a)/n)^(n+1))/(4*(n+1))
%%
f = @(x) sqrt(1+x);
I = trapcomp(-1,3,4,f)
%%
e1 = 1e-1;
M1 = 10;
M2 = 100;
e2 = e1*(M1/M2)^4
%%
f = @(x,y) exp(x+3*y);
a= 0; b = 1; c = 0; d = 1;
eps = [-1/sqrt(3) 1/sqrt(3)];
x = (a+b)/2 + (b-a)/2 .* eps;
y = (c+d)/2 + (d-c)/2 .* eps;
S = 0;
for i = 1:2
    for j = 1:2
        S = S + f(x(i),y(j));
    end
end
I = (b-a)*(d-c)*S/4
%%
f = @(t,y) -(1+t)*y;
y0 = 2;
Heun(f,10,2)
%%
A = @(y) diagonals([3 -2 -1],9);
g = @(t) 0;
fun = @(t,y) A(y)*y + g(t);
h = linspace(0,1,1000);
l = max(eig(A(1)));
z = h.*l;
R = 1 + z;
R_abs = real(R).^2 + imag(R).^2;
plot(h,R_abs)
yline(1)
h = 2*real(l)/(real(l).^2+imag(l).^2);
%%
A = @(y) diagonals([3 -2 -1],9);
g = @(t) exp(sin(pi.*t));
fun = @(t,y) A(y)*y + g(t);
y0 = 4*ones(9,1);
[x,u] = eulero_avanti_sistemi(fun,[0 10],y0,10/0.1);
u(5,2)
u(5,end)
[u5_min ,i] = min(u(5,:));
u5_min
x(i)
h = [1e-2 5*1e-3 2.5*1e-3 1.25*1e-3];
err = [];
y_es = 1.142174435142178;
for i = h
    [x,u] = eulero_avanti_sistemi(fun,[0 10],y0,10/i);
    err = [err max(abs(y_es-u(5,end)))];
end
p = stimap_2(err,h);
p(end)
%%
A = @(y) diagonals([3 -2 -1],9);
g = @(t) exp(sin(pi.*t));
fun = @(t,y) A(y)*y + g(t);
y0 = 4*ones(9,1);
[x, u] = Crank_Nicolson_sistemi(fun,[0 10],y0,10/0.1);
plot(x,u)
%%
A = @(y) diagonals([3 -2 -1],9);
g = @(t) exp(sin(pi.*t));
fun = @(t,y) A(y)*y + g(t);
y0 = 4*ones(9,1);
A_RK = [1/3 0; 1/3 1/3];
b = [1/2 1/2];
c = [1/3 2/3];
[x,u] = Runge_stronzo(A_RK,b,c,fun,[0 10],y0,10/0.1);
%%
x = 0:4;
y = [3 1 3 4 5];
polyval(polyfit(x,y,4),1.5)
%%
f = @(x) exp(x);
poly_stima(f,4,0,2,'equi')
%%
f = @(x) x.^3 + sin(pi.*x);
[y,x,y0] = interp_tratti(f,0,3,2,3,'equi',1.5);
plot(x,f(x),x,y)
y0
%%
f = @(x) x.^2 + abs(x);
x_nod = -1:2/4:1;
x = -1:0.0001:1;
y = polyval(polyfit(x_nod,f(x_nod),2),x);
plot(x,f(x),x,y,'LineWidth',2)
trapz(x,y)
%%
e1 = 1e-1;
M1 = 10;
M2 = 100;
e2 = e1*(M1/M2)^2
%%
syms t y;
f = -t^2*y;
EAI_sym(f,1,1)
%%
f = @(t,y) -(t+1).*y.^3;
y0 = 1;
A = [0 0; 1/2 0];
b = [0 1];
c = [0; 1/2];
[t,u] = Runge_Kutta_gen(A,b,c,f,10,0.1,1);
plot(t,u)
u(end)
%%
A = [-2 1; 2 -3];
y0 = [1 3]';
fun = @(t,y) A*y;
[t,u] = Crank_Nicolson_sistemi(fun,[0 1],y0,10);
u(:,end)
%%
A = @(y) [-4 -100; 1 0];
g = @(t) [200*cos(10*t); 0];
y0 = [50 0]';
tv = [0 5];
fun = @(t,y) A(y)*y + g(t);
[t,u] = eulero_avanti_sistemi(fun,tv,y0,5/1e-2);
plot(t,u(2,:),'LineWidth',2);
[u_min,i] = min(u(2,:));
x = @(t) 5.*sin(10.*t);
h = [1e-3 5*1e-4 2.5*1e-4 1.25*1e-4];
err = [];
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,5/i);
    err = [err max(abs(x(t(end))-u(2,end)))];
end
err
p = stimap_2(err,h);
p(end)
A = [-4 -100; 1 0];
l = max(eig(A));
h = 0:0.0001:1;
z = h.*l;
R = 1 + z;
R_abs = real(R).^2+imag(R).^2;
plot(h,R_abs)
yline(1)
ylim([0 2])
f = @(h) real(1+h.*l).^2+imag(1+h.*l).^2-1;
fsolve(f,0.5)
%%
K = 1;
A = @(y) [-1 -(147/16+K*y(2)); 1 0];
g = @(t) [-9*exp(-t/2)*(sin(3*t)^2-1)-9/2*exp(-t/4)*sin(3*t),0]';
tf = 5;

fun = @(t,y) A(y)*y + g(t);
tv = [0 tf];
y0 = [-3/4 3]';

h = 1e-2;
Nh = tf/h;
[t, u] = eulero_avanti_sistemi(fun, tv, y0, Nh);
plot(t,u,'LineWidth',2);
u(:,2)
u(:,end)
[u_min, i] = min(u(2,:));
u_min
t(i)

y_es = @(t) (3*exp(-t/4)).*[-1/4*cos(3*t)-3*sin(3*t); cos(3*t)];
subplot(2,1,1);
plot(t,u,'LineWidth',2)
subplot(2,1,2);
plot(t,y_es(t),'LineWidth',2)

h = [1e-3 5*1e-4 2.5*1e-4 1.25*1e-4];
err = [];
for i = h
    Nh = tf/i;
    [t, u] = eulero_avanti_sistemi(fun, tv, y0, Nh);
    err = [err norm(u(:,end)-y_es(t(end)))];
end
err;
p = log(err(end)/err(end-1))/log(h(end)/h(end-1));
close all
h = 0:0.001:1;
A = [-1 -(147/16); 1 0];
l = max(eig(A));
z = @(h) h.*l;
R = @(h) 1 + z(h) + z(h).^2/2;
R_abs = @(h) (real(R(h)).^2+imag(R(h)).^2 - 1);
fsolve(R_abs,0.5)
%%

% dy/dt = A*y
% y = A\f(t,y)
K = 0;
A = @(y) [-1 -(147/16+K*y(2)); 1 0];
g = @(t) 0;
tf = 5;
fun = @(t,y) A(y)*y + g(t);
tv = [0 tf];
y0 = [-3/4 3]';
[t,u] = Crank_Nicolson_sistemi(fun,tv,y0,5/1e-2);
plot(t,u,'LineWidth',2)
u(:,2)
u(:,end)
%%
f = @(x) sin(x + sqrt(2));
x_nod = 0:(pi/3):pi;
P = polyfit(x_nod,f(x_nod),3);
dP = polyder(P);
polyval(P,2)
polyval(dP,2)
%%
f = @(x) exp(3.*x);
st = poly_stima(f,3,-1,1,'CGL')
%%
f = @(x) ((x./pi).^2).*(x>= 0 & x< pi) + ((2-pi./x).^2).*(x>= pi & x < 2*pi);
x_nod = 2*pi/5.*(0:5);
cubicspline(x_nod,f(x_nod),x_nod(4))
%%
x = 0:4;
y = [2 2 0 1 0];
poly_minim_IMQ(1,x,y)
%%
f = @(x) exp(x);
a = -1;
b  = 3;
 I = trapcomp(a, b, 4, f)
%%
f = @(x) x.^3;
F = @(c) 3.*(f(3+c)+f(3-c))-6.^4/4;
fsolve(F,1)
%%
f = @(t,y) -2*y;
y0 = 7;
[t,u] = teta_met(f,[0 5],y0,5/(1e-1),3/4);
plot(t,u,'LineWidth',2)
%%
tf = 10;
y0 = [2 0]';
[t,u] = eulero_avanti_sistemi(@prova,[0 tf],y0,tf/1e-2);
u(2,2)
u(2,end)
figure(1)
plot(t,abs(u(2,:)))
yline(0.1)
figure(2)
subplot(2,1,1)
plot(t,u(2,:))
x = @(t) 2*exp(-t./2).*sin(t);
subplot(2,1,2)
plot(t,x(t))

h = [1e-3 5*1e-4 2.5*1e-4 1.25*1e-4];
err = [];
for i = h
    [t,u] = eulero_avanti_sistemi(@prova,[0 tf],y0,tf/i);
    err = [err max(abs(x(t)-u(2,:)))];
end
err
p = stimap_2(err,h);
2/1.04

t = 0:0.001:10;
df =@(x) [-2 1; -20*(x) 0];
for i = 1 : length(t)
    l(i) = max(eig(df(x(t(i)))));
end
l = max(l)
h = 0:0.001:1;
R = @(h) 1 + h*l;
sol = @(h) real(R(h)).^2 + imag(R(h)).^2-1;
fsolve(sol,0.1)
close all
%%
tf = 10;
y0 = [2 0]';
tv = [0 tf];
[t,u] = multipasso(@prova,tv,y0,10/1e-2);
u(2,2)
u(2,3)
u(2,end)
h = [1e-3 5*1e-4 2.5*1e-4 1.25*1e-4];
err = [];
for i = h
    [t,u] = multipasso(@prova,[0 tf],y0,tf/i);
    err = [err max(abs(x(t)-u(2,:)))];
end
err
p = stimap_2(err,h);
%%
x = linspace(0,1,50);
rng(1);
y = 2 * x.^2 + 0.2 * sin(100*pi*x) + 0.2 * randn(1,50);
y_IMQ = polyval(polyfit(x,y,2),0.5)
%%
h = 10/5;
x_nod = 0:h:10;
f = @(x) (x+2).^2;
interp1(x_nod,f(x_nod),1)
%%
f = @(x) sin(10.*x)-7.*x;
df = Jac(f);
ddf = Jac(df);
a = 0;
b = 3;
tol = 1e-4;
x = 0:0.001:b;
H = ceil(3/sqrt(tol*8/max(ddf(x))))
%%
n = 6;
a = -1;
b = 1;
str = 'CGL';
[eq, phi, A] = poly_equation(n,a,b,str);
w = matlabFunction(eq);
pmedcomp(a,b,1,w);
x = a:0.01:b;
plot(x,w(x))
%%
f = @(x) exp(x);
integr_stima(-2,2,1e-1,f,'t')
%%
f = @(t,y) -2*(t+1)*y^2;
y0 = 3;
A = [0 0; 3/4 0];
b = [1/3 2/3];
c = [0 3/4];
[t,u] = Runge_Kutta_gen(A,b,c,f,10,0.1,y0);
plot(t,u)
%%
f = @(t,y) -5*(t+1)*y;
y0 = 8;
tf = 100;
[t,u] = Crank_Nicolson(f,1,y0,0.2);
df = @(t) -5*(t+1);
t1 = 0:0.001:100;
l = max(df(t1));
h = 0:0.001:1;
z = h.*l;
R = (1+z./2)./(1-z./2);
R_abs = real(R).^2 + imag(R).^2;
figure(1)
plot(h,R_abs)
ylim([0 2])
yline(1)
figure(2)
plot(t,u)
%%
A = [-3 -4; 1 0];
fn = @(t,y) A*y;
y0 = [7;1];
l = max(eig(A));
h = 0:0.001:1;
R_abs = @(h) abs(1 + h.*l).^2 - 1;
fsolve(R_abs,0.5)
%%
A = @(y) [-2 -6; 1 0];
g = @(t) [10 0]';
y0 = [1 4]';
tv = [0 5];
fun = @(t,y) A(y)*y + g(t);
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tv(end)/0.1);
u(2,2);
u(2,end);
x = @(t) exp(-t).*(7/3.*cos(sqrt(5).*t) + 10/(3*sqrt(5)).*sin(sqrt(5).*t))+5/3;
figure(1)
subplot(2,1,1)
plot(t,u(2,:),'LineWidth',2)
subplot(2,1,2)
plot(t,x(t),'LineWidth',2)
h = [1e-2 5*1e-3 2.5*1e-3 1.25*1e-3];
err = [];
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tv(end)/i);
    err = [err abs(x(t(end))-u(2,end))];
end
figure(2)
loglog(h,err,'r',h,h,'--k',h,h.^2,'-.k','LineWidth',2);
legend('Eh','H','H^2');
stimap_2(err,h);

[t,u] = Heun_sistemi(fun,tv,y0,tv(end)/0.1);
figure(3)
plot(t,u(2,:),'LineWidth',2)
u(2,2);
u(2,end);
close all
l = max(eig(A(1)));
h = 0:0.001:1;
R = abs(1 + h.*l + (h.*l).^2/2).^2;
plot(h,R)
yline(1)
ylim([0 2])
R_fun = @(h) abs(1 + h.*l + (h.*l).^2/2).^2 -1;
fsolve(R_fun,1)
%%
f = @(x) cos(3.*x + sqrt(5));
poly_stima(f,4,0,1,'equi')
%%
f = @(x) exp(x);
[y,x,y0] = interp_tratti(f,-2,2,2,4/1,'equi',1.5);
y0
plot(x,f(x),x,y)
%%
f = @(x) exp(x./10) + 0.1*sin(pi.*x + sqrt(2));
x_nod = 0:10;
polyval(polyfit(x_nod,f(x_nod),2),11)
%%
f = @(x) exp(x);
a = [5/9 8/9 5/9];
x_nod = [-sqrt(3/5) 0 sqrt(3/5)];
I = sum(a.*f(x_nod))
gausscomp_2(-2,2,1,f)
%%
f = @(t,y) 3*exp(y) -19*t;
y0 = 0;
[t,u] = eulero_avanti(f,0.1,0,0.1)
%%
h = 10/5;
x_nod = 0:h:10;
f = @(x) x.^3 + 3;
x = 0:0.01:10;
interp1(x_nod,f(x_nod),3)
%%
e1 = 1e-2;
M1 = 10;
e2 = 1e-5;
M2 = ceil(M1*nthroot(e1/e2,4))
%%
f = @(x) sin(pi*x)+2*x^4-7*x^3;
integr_stima(-2,2,1e-4,f,'s')
%%
syms b;
dddf = 6*b;
h = 1/4;
expand(-1/12*h^2*2*dddf)
%%
f = @(t,y) -exp(y) + 2*t;
h = 0.1;
y0 = 3;
[t,u] = Heun(f,10,y0,0.1);
plot(t,u)
%%
A = @(y) [-2 -8 0; 1 0 -1; 0 0 -1];
g = @(t) exp(-t/2).*[(-3/2*(4*pi*sin(pi*t)+(4*pi^2-29)*cos(pi*t))) 2 1]';
y0 = [-3 6 2]';
tf = 10;
fun = @(t,y) A(y)*y + g(t);

tv = [0 tf];
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/0.1);
u(2,2)
u(2,end)
plot(t,u,'LineWidth',2)
l = max(eig(A(1)));
h = 0:0.001:1;
R_abs = @(h) abs(1 + h*l)-1;
[t,u] = Crank_Nicolson_sistemi(fun,tv,y0,tf/0.1);
u(2,2)
u(2,end)
y = @(t) exp(-t/2).*[-3*cos(pi*t)-6*pi*sin(pi*t),6*cos(pi*t),2]';
h = [1e-2 5*1e-3 2.5*1e-3 1.25*1e-3];
err = [];
for i = h
    [t,u] = Crank_Nicolson_sistemi(fun,tv,y0,tf/i);
    err = [err norm(u(:,end)-y(tf))];
end
err
p = stimap_2(err,h);
loglog(h,err,h,h,'--k',h,h.^2,'-.k','LineWidth',2)
legend('Eh','H','H^2')
%% 
f = @(x) 5 + x.^7;
gausslegendre_comp(0,2,f,1,2,'equi')
%%
f = @(x) abs(1-exp(abs(sin(x))));
x = -2:0.001:2;
y = poly_lagr(f,4,-2,2,x,'equi');
max(y)
%%
f = @(x) x.*exp(sin(x));
h = 10/10;
x_nod = 0:h:10;
x_dis = 0:0.01:10;
y = interp1(x_nod,f(x_nod),x_dis);
[e_max,i] = max(abs(y-f(x_dis)));
e_max
x_dis(i)
%%
b = 1;
a = -1;
ddf = exp(1);
toll = 1e-2;
N = ceil(sqrt(((b-a)^3)*ddf/(24*toll)))
%%
f = @(x) 5 + x.^7;
I = gausslegendre_comp(0,3,f,3,1,'equi')
%%
syms h;
e = simplify(-84/12*h)
%%
f = @(t,y) -(1+t)*y;
y0 = 3;
[t,u] = Crank_Nicolson(f,10,3,0.1);
plot(t,u)
u(2)
%%
a = -1/4; b = -1/5; c = -1/8; d = 0;
A = @(y) [a b*y(1); c*y(2) d];

g1 = @(t) 8*exp(-t/4)*(cos(pi*t)^2 - pi*sin(pi*t));
g2 = @(t) 5*(exp(-t/4)*cos(pi*t)^2 - pi*sin(pi*t));
g = @(t) [g1(t); g2(t)];

y1 = 8; y2 = 5;
y0 = [y1 y2]';

fun = @(t,y) A(y)*y + g(t);
tf  = 10;
tv = [0 tf];
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/0.05);
plot(t,u(2,:),'LineWidth',2)
u(:,end);
[t,u] = metodo_numerico(A,g,tv,y0,tf/0.05);
u(:,end);
plot(t,u);
y = @(t) [8*exp(-t/4)*cos(pi*t),5*cos(pi*t)]';

h = [1e-3 5*1e-4 2.5*1e-4 1.25*1e-4];
err = [];
for i = h
    [t,u] = metodo_numerico(A,g,tv,y0,tf/i);
    err = [err norm(u(:,end)-y(t(end)))];
end
err
stimap_2(err,h);
loglog(err,h,h,h,'--k',h,h.^2,':k','LineWidth',2);
legend('Eh','H','H^2')
%%
f = @(x) x + sin(pi*x + 1);
h = 5/10;
x_nod = 0:h:5;
x_dis = 0:0.001:6;
y_dis = polyval(polyfit(x_nod,f(x_nod),1),x_dis);
plot(x_nod,f(x_nod),'o',x_dis,y_dis,'LineWidth',2)
polyval(polyfit(x_nod,f(x_nod),1),6)
%%
f = @(x) exp(x)
a = 0;
b = 2;
tol = 1e-2;
H = sqrt(8*tol/f(2));
M = ceil(2/H)
x_nod = a:2/M:b;
x_dis = a:0.001:b;
err = max(abs(f(x_dis)-interp1(x_nod,f(x_nod),x_dis)));
%%
x_nod = 1:5;
y_nod = [1 1 0 1 2];
x_dis = 1:0.001:5;
y_dis = spline(x_nod,y_nod,x_dis);
plot(x_nod,y_nod,'o',x_dis,y_dis);
spline(x_nod,y_nod,4.5)
%%
syms b g;
f = @(x) b.*x.^3 + g.*x.^2 + 1;
I = pmedcomp(-2,2,2,f)
%%
f = @(x) 2.^x;
x = 0;
h = 1/4;
df = (-3*f(x)+4*f(x+h)-f(x+2*h))/(2*h)
%%
f = @(t,y) 1 + t*y;
y0 = 1;
[t,u] = Crank_Nicolson(f,2,1,0.1);
plot(t,u)
u(2)
%%
c = 2;
k = @(t) 3*(2-exp(-t));
z = @(t) 10*exp(-t)*cos(t)*(30*exp(-t)*(2-exp(-t))*cos(t)-2);
m = 1;
fun = @(t,y) [-c/m*y(1)-k(t)/m*y(2)^2 + z(t); y(1)];
v0 = -10; x0 = 10; tf = 2;
tv = [0 tf];
y0 = [v0; x0];

[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/0.05);
figure(1)
plot(t,u,'LineWidth',2)
u(:,end);

A = @(t,y) [-c/m -k(t)/m*y(2); 1 0];
g = @(t) [z(t);0];
[t,u] = metodo_numerico(A,g,tv,y0,tf/0.05);
figure(2)
plot(t,u,'LineWidth',2)
u(:,end);

y = @(t) [-10*exp(-t)*(cos(t)+sin(t));10*exp(-t)*cos(t)];
h = [10 5 2.5 1.25].*1e-4;
err = [];
for i = h
    [t,u] = metodo_numerico(A,g,tv,y0,tf/i);
    err = [err norm(u(:,end)-y(t(end)))];
end
format("shortE")
err
stimap_2(err,h)
loglog(h,err,h,h,'--k',h,h.^2,':k');
legend('Eh','H','H^2')
%%
a = 0; b = 1;
alpha = 1; beta = 0;
mu = 1;
sigma = @(x) 5 + 0.*x;
f = @(x) 1 + 0.*x;
N = 10-1;
h = 1/10;
x_nodes = linspace(a,b,N+2);

A = sparse(1:N,1:N,2,N,N) + ...
    sparse(2:N,1:N-1,-1,N,N) + ...
    sparse(1:N-1,2:N,-1,N,N);
A = mu/h^2 * A;
A = A + sparse(diag(sigma(x_nodes(2:end-1))));
bv = (f(x_nodes(2:end-1))');
bv(1) = bv(1) + alpha*mu/h^2;
bv(end) = bv(end) + beta*mu/h^2;

u_h = A\bv;
u_h = [alpha; u_h; beta];
plot(x_nodes,u_h)
u_h(6)
%%
A = @(y) diagonals([3 -2 -1],9);
g = @(t) exp(sin(pi.*t));
y0 = 4*ones(9,1);
tf = 10;
tv = [0 tf];
fun = @(t,y) A(y)*y + g(t);
[t,u1] = eulero_avanti_sistemi(fun,tv,y0,tf/0.1);
subplot(2,1,1)
plot(t,u1,'LineWidth',2)
A_RK = [0];
b = [1];
c = [0];
[t, u2] = Runge_stronzo(A_RK,b,c,fun, tv, y0, tf/0.1);
subplot(2,1,2)
plot(t,u2,'LineWidth',2)

if u1 == u2
    disp('vettori uguali con EA')
end
[t,u1] = Heun_sistemi(fun,tv,y0,tf/0.1);
A_RK = [0 0; 1 0];
b = [1/2 1/2];
c = [0 1]; 
subplot(2,1,1)
plot(t,u1,'LineWidth',2)
[t, u2] = Runge_stronzo(A_RK,b,c,fun, tv, y0, tf/0.1);
plot(t,u2,'LineWidth',2)
if u1 == u2
    disp('vettori uguali con Heun')
end

%%
fun = @(t,y) -(t+1)*y^3;
y0 = 1;
tv = [0 10];
A = [0 0; 1/2 0];
b = [0 1];
c = [0 1/2];
[t1,u1] = Runge_stronzo(A,b,c,fun,tv,y0,10/0.1);
[t2,u2] = eulero_avanti_sistemi(fun,tv,y0,tf/0.1);
plot(t1,u1,t2,u2,'LineWidth',2)
legend('RK','EA');
u(2)
%%
K = 1;
A = @(y) [-1 -(147/16+K*y(2)); 1 0];
g = @(t) [-9*exp(-t/2)*(sin(3*t)^2-1)-9/2*exp(-t/4)*sin(3*t);0];
fun = @(t,y) A(y)*y + g(t);
tf = 5;
tv = [0 tf];
y0 = [-3/4 3]';
subplot(2,1,1)
[t,u1] = eulero_avanti_sistemi(fun,tv,y0,tf/1e-2);
u1(:,end)
plot(t,u1,'LineWidth',2)
subplot(2,1,2)
A = [0 0; 2/3 0];
b = [1/4 3/4];
c = [0 2/3];
[t,u2] = Runge_stronzo(A,b,c,fun,tv,y0,tf/1e-2);
plot(t,u2,'LineWidth',2)
u2(:,end)
%%
fun = @(t,y) -2.*(1+sin(pi.*t)).*(2-y).^2;
tf = 5;
tv = [0 tf];
y0 = 3;
h = 1/4;
A = [0 0 0 0;...
    1/2 0 0 0;...
    0 1/2 0 0;...
    0 0 1 0];
c = [0 1/2 1/2 1];
b = [1/6 1/3 1/3 1/6];

[t1,u1] = Runge_stronzo(A,b,c,fun,tv,y0,tf/h);
[t2,u2] = Runge_Kutta_4(fun,tf,y0,h);
subplot(2,1,1)
plot(t1,u1,'LineWidth',2)
subplot(2,1,2)
plot(t2,u2,'LineWidth',2)

%%
K = 1;
A = @(y) [-1 -(147/16+K*y(2)); 1 0];
g = @(t) [-9*exp(-t/2)*(sin(3*t)^2-1)-9/2*exp(-t/4)*sin(3*t);0];
fun = @(t,y) A(y)*y + g(t);
tf = 5;
tv = [0 tf];
y0 = [-3/4 3]';
A = [1 1; 1 1];
b = [1/4 3/4];
c = [0 2/3];
[t,u2] = Runge_Kutta(A,b,c,fun,tv,y0,tf/1e-2);
plot(t,u2,'LineWidth',2)
%%
z = @(t) exp(-t/2)*(2*cos(t)-7/2*sin(t))+40*exp(-t)*sin(t)^2;
f = @(y) [y(2); -10*y(1)^2-2*y(2)];
g = @(t) [0;z(t)];
fun = @(t,y) f(y)+g(t);
tf = 10;
tv = [0 tf];
y0 = [0;2];
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/1e-2);
plot(t,abs(u(1,:)),'-o','LineWidth',1)
yline(0.1)
i = find(abs(u(1,:))>=0.1,1,'last');
t_m = t(i+1);

x_es = @(t) 2.*exp(-t./2).*sin(t);
h = [10 5 2.5 1.25].*1e-4;
err = [];
for i = h
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    err = [err max(abs(u(1,:)-x_es(t)))];
end
err
stimap_2(err,h);

A = @(x) [0 1; -20*x -2];
t = 0:0.001:10;
for i = 1:length(t)
    l(i) = max(eig(A(x_es(t(i)))));
end
l = max(l)
h = 0:0.001:1;
z = h.*l;
R_abs = abs(1 + z);
plot(h,R_abs,'LineWidth',2);
yline(1)
ylim([0 2])
R_abs = @(h) abs(1+h*l)-1;
fsolve(R_abs,1)
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/0.01);
plot(t,u(1,:))
%%
f = @(x) exp(x./2);
a = 0;
b = 3;
h = (b-a)/3;
x_nodes = a:h:b
n = 2;
h2 = (1-0)/n;
x_nodes2 = 0:h2:1;
polyval(polyfit(x_nodes2,f(x_nodes2),n),0.75)

%%
a = 0; b = 1;
alpha = 0; beta = 0;
Nh = 2-1;
x_nodes = linspace(a,b,Nh+2);
h = 1/2;
N = Nh;
mu = @(x) 1.*(x>0 & x < 1/2) + 2.*(x>=1/2 & x<1);
fun = @(u) -1/h^2*(mu(0.5+1/4)*(0-u)-mu(0.5-1/4)*(u-0))-5;
fsolve(fun,1)
%%
n = 5; a = -1; b = 1;
x = a:0.001:b;
i = 0:n;
x_nod = -cos(pi.*i./n);
syms y;
omega = 1;
for i = 1 : n + 1
    omega = omega*(y-x_nod(i));
end
omega = matlabFunction(simplify(expand(omega)));
stim = abs(1e3)*max(abs(omega(x)))/factorial(n+1)
%%
A = diagonals([4 -2 -2],10);
g = @(t) cos(pi.*t);
fun = @(t,y) A*y + g(t);
y0 = 7*ones(10,1);
tf = 10;
tv = [0 tf];
h = 0.1;
subplot(2,1,1)
[t,u] = eulero_indietro_sistemi(fun,tv,y0,tf/h);
plot(t,u)
u(3,2)
u(3,3)

subplot(2,1,2)
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
plot(t,u)
u(3,2)
u(3,end)

