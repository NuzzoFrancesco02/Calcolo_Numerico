%%
x = [0 0.5 1];
y = [2 1 1.5];
P = polyfit(x,y,2);
dP = polyder(P);
polyval(P,0.25)
polyval(dP,0.25)
%%
x = 0:0.25:1;
y = [3 0.5 1.5 -0.5 1];
x_dis = 0:0.001:1;
y = interp1(x,y,x_dis);
min(y)
%%
f = @(x) sin(pi.*x);
n = 4;
h = 2/4;
x_nod = 0:h:2;
interp1(x_nod,f(x_nod),1.6)
%%
x = [0:0.25:1];
y = [3 0.5 1.5 -0.5 1];
poly_scarto_IMQ(x,y,2)
%%
f = @(x) sqrt(2+abs(x));
I = trapcomp(-1,2,3,f)
%%
M1 = 20;
M2 = ceil(M1*sqrt(1e-1/1e-3))
%%
f = @(t,y) -(1+sin(t)).*(y.^2)./81;
df = @(t,y) -2*(1+sin(t)).*(y)./81;
y0 = 9;
h = 0.2;
[t,u] = eulero_indietro_newton(f,df,5,y0,h);
u(end)
%%
h = 1e-2;
Nh = 5/h;
[t,u] = eulero_avanti_sistemi(@prova,[0 5],[-3/4 3]',Nh);
u(:,2)
u(:,end)
[u_min,i] = min(u(2,:))
u_min
t(i)
y = @(t) 3*exp(-t/4).*[(-cos(3*t)/4-3*sin(3*t)),cos(3*t)]';
h = [1e-3 5*1e-4 2.5*1e-4 1.25*1e-4];
err = [];
for i = h
    Nh = 5/i;
    [t,u] = eulero_avanti_sistemi(@prova,[0 5],[-3/4 3]',Nh);
    err=[err norm(y(t(end))-u(:,end))];
end
err
p = stimap_2(err,h);
p(end)
A = [-1 -(147/16);1 0];
t = 0:0.001:5;

2/max(abs(eig(A)))
%%
f = @(x) sin(x+sqrt(2));
n = 3;
h = pi/n;
x_nod = 0:h:pi;
P = polyfit(x_nod,f(x_nod),3);
dP = polyder(P);
polyval(P,2)
polyval(dP,2)
%%
f = @(x) exp(3*x);
poly_stima(f,3,-1,1,'CGL')
%%
f = @(x) ((x./pi).^2).*(x>=0 & x < pi)+( ...
   (2-x./pi).^2).*(x>=pi & x <2*pi);
n = 5;
i = 0:n;
x_nod = 2*pi.*i./n;
cubicspline(x_nod,f(x_nod),x_nod(4))
%%
x = 0:4;
y = [2 2 0 1 0];
poly_minim_IMQ(1,x,y)
%%
f = @(x) exp(x)
a = -1; b = 3;
H = 1;
N = (b-a)/H;
I = trapcomp(a,b,N,f)
%%
x = 6;
f = @(x) x.^3;
fun = @(c) 3*(f(3+c)+f(3-c))-(x^4)/4;
x = -1:0.01:3;
fsolve(fun,2)
%%
f = @(t,y) -2*y;
y0 = 7;
teta = 3/4;
[t,u] = teta_met(f,[0,0.1],y0,1,teta);
u(2)
%%
h = 1e-2;
Nh = 10/h;
[t,u] = eulero_avanti_sistemi(@prova,[0 10],[2 0]',Nh);
u(2,2)
u(end,end)
plot(t,u(2,:),t,u(1,:))
hold on
yline(0.1)

y = @(t) 2.*exp(-t./2).*sin(t);
h = [1e-3 5*1e-4 2.5*1e-4 1.25*1e-4];
err = [];
for i = h
    Nh = 10/i;
    [t,u] = eulero_avanti_sistemi(@prova,[0 10],[2 0]',Nh);
    err = [err max(abs(u(2,:)-y(t)))];
end
err
%%
z = @(t) exp(-t./2).*(2.*cos(t)-7.*sin(t)./2)+40.*exp(-t).*sin(t).^2;
fn = @(t,y1,y2) [z(t)-2*y1-10*y2^2,y1]';
y = @(t) 2.*exp(-t./2).*sin(t);
t = 0:0.001:100;
for i = 1 : length(t)
    J = [-2 1;-40*exp(-t(i)./2).*sin(t(i)) 0];
    l(i) = max(abs(eig(J)));
end
l = max(l);
2/l
%%
f = @(x) sin(x+sqrt(5));
n = 3;
h = (pi/n);
x = 0:0.001:pi;
[y, P] = poly_lagr(f,3,0,pi,x,'CGL');
dP = polyder(P);
polyval(P,1)
polyval(dP,1)
plot(x,y,x,f(x))
%%
f = @(x) sin(pi.*x)
poly_stima(f,3,-1,1,'equi')
%%
0
%%
x = 0:4;
y = [1 1 2 2 3];
polyval(polyfit(x,y,1),4)
%%
f = @(x) exp(x);
I = pmedcomp(-2,3,5,f)
%%
e1 = 0.08;
e2 = e1*(25/50)^2
%%
f = @(t,y) -3*y;
y0 = 4;
[t,u] = teta_met(f,[0 0.1],y0,1,1);
u(end)
%%
f = @(x) exp(x).*cos(x);
a = -2;
b = 2;
h = [1 0.5 0.1 0.05 0.01];
x = a:0.001:b;
err = [];
for i = h
    x_nod = a:i:b;
    y = cubicspline(x_nod,f(x_nod),x);
    err = [err max(abs(f(x)-y))];
end
stimap_2(err,h)
%%
f = @(x) exp(2.*x + sin(pi.*x));
x = -1:0.0001:1;
y = poly_lagr(f,5,-1,1,x,'equi');
[err_max, i] = max(abs(f(x)-y));
err_max
x(i)
%%
x = 0:0.25:1;
y = [1 0.5 1.5 -0.25 1];
x_dis = 0:0.001:1;
y_dis = interp1(x,y,x_dis);
plot(x_dis,y_dis)
%%
f = @(x) 2*abs(sin(pi.*x));
[y,x,y0] = interp_tratti(f,0,4,2,4,'equi',1.75);
plot(x,y,x,f(x))
y0
%%
x = [0 1/2 2];
phi = poly_phi(x);
phi0 = matlabFunction(phi(1));
simpcomp(0,2,1e5,phi0)
%%
E = [1e-1 1e-4];
H = [10 100]
i = 2;
log(E(i)/E(i-1))/log(H(i)/H(i-1))
%%
f = @(x,y) 2.^(x+3*y);
a = 0; b = 1; c = 0; d = 1;
I = (b-a)*(d-c)*(f(a,c)+f(a,d)+f(b,c)+f(b,d))/4
%%
f = @(t,y) -sqrt(10*y/(10+t));
y0 = 9;
[t,u ]= Heun(f,4,y0,0.2);
plot(t,u)
u(end)
%%
f = @(x) 1 +x.^5 + abs(x);
trapcomp(-1,1,1,f)
%%
syms t h y;
f = -t*y;
y0 = 5;
EAI_sym(f,1,y0)
%%
f = @(t,y) -(t+1)*y.^2;
t0 = 0;
y0 = 9;
A = [0 0; 1 0];
c = [0 1]';
b = [1/2 1/2];
u = Runge_Kutta_gen(A,b,c,f,0.1,0.1,9)
%%
f = @(t,y,v) -2*y;
y0 = 10;
w0 = 0;
h = 0.1;
[t,u,v] = Leap_Frog(f,[0 5],h,y0,w0);
plot(t,u,t,v)  
%%
[t,u] = eulero_avanti_sistemi(@prova,[0 100],[-5/12 10]',100/(1e-2));
u(2,2)
u(2,3)
u(2,end)
%plot(t,u(2,:))
f = @(t) 10*exp(-t./24).*cos(sqrt(2303/576).*t);
plot(t,f(t))