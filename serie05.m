%% serie 5
f = @(x) x.*sin(x);
%% 1.1
a = -2;
b = 6;
x_nod = linspace(a,b,1000);
P_dis =f(x_nod);
figure()
plot(x_nod,P_dis,'LineWidth',2);
%% 1.2
PP_dis = [];
err_dis = [];
err_inf = [];
for n = 2:2:6
    h = (b-a)/n;
    x_nod = a:h:b;
    f_nod = f(x_nod);
    P = polyfit(x_nod,f_nod,n);
    P_dis = polyval(P,x_nod);
    PP_dis = [PP_dis;P_dis];
    err_dis = [err_dis;abs(P_dis-P_dis)];
    err_inf = [err_inf; max(err_abs)];
end
figure()
plot(x_nod,P_dis,'LineWidth',2)
hold on
plot(x_nod,PP_dis,'LineWidth',2)
figure()
plot(x_nod,err_dis,'LineWidth',2)
%% 2.1 2.2
clear 
clc
f = @(x) sin(1./(1+x.^2));
a = -2*pi;
b = 2*pi;
n = [2 4 8 10];

x_nod = linspace(a,b,1000);
P_dis = f(x_nod);
PP_dis = [];
err_inf = [];
for n = n
    h = (b-a)/n;
    x_nod = a:h:b;
    P = polyfit(x_nod,f(x_nod),n);
    P_dis = polyval(P,x_nod);
    PP_dis = [PP_dis;P_dis];
    err_inf = [err_inf;abs(P_dis-P_dis)];
end
figure()
plot(x_nod,P_dis,'LineWidth',2);
hold on;
plot(x_nod,PP_dis,'LineWidth',2);
figure()
plot(x_nod,err_inf,'LineWidth',2);

%% 2.3
clear 
clc
f = @(x) sin(1./(1+x.^2));
a = -2*pi;
b = 2*pi;
n = [2 4 8 10];

x_nod = linspace(a,b,1000);
P_dis = f(x_nod);
PP_dis = [];
err_inf = [];
for n = n
    k = 0 : n;
    t = -cos(pi.*k./n);
    x_nod = ((b-a).*t./2) + (a+b)./2;
    P = polyfit(x_nod,f(x_nod),n);
    P_dis = polyval(P,x_nod);
    PP_dis = [PP_dis;P_dis];
    err_inf = [err_inf;abs(P_dis-P_dis)];
end
figure()
plot(x_nod,P_dis,'LineWidth',2);
hold on;
plot(x_nod,PP_dis,'LineWidth',2);
figure()
plot(x_nod,err_inf,'LineWidth',2);

%% 3
clear
clc
f = @(x) -x.^3 +3.*x.^2 - 2;
x_dis = [0 0.5 2];
n = length(x_dis)-1;
phi = [];
for k = 1 : n
    i = 0:n;
    xi = x_dis(i+1);
    xk = x_dis(k);
    phi_f = @(x) prod((x-xi)./(xk-xi));
 
end
Prod = product(phi);
P_n = @(x) sum(f(x).*Prod(x));

plot(x_dis,f(x_dis));
hold on
plot(x_dis,P_n(x_dis))

%%
clear
clc
sigma = [0.18 0.3 0.5 0.6 0.72 0.75 0.8 0.9 1];
epsilon = [0.0005 0.001 0.0013 0.0015 0.0020 0.0045 0.006 0.007 0.0085];
n = length(epsilon)-1;
figure()

x_dis = linspace(min(sigma),max(sigma),1000);

subplot(3,2,1)
P = polyfit(sigma,epsilon,n);
y_lag = polyval(P,x_dis);
plot(sigma,epsilon,'o','LineWidth',2);
hold on
plot(x_dis,y_lag,'LineWidth',2)
title('Lagrange')

subplot(3,2,2)
y_interp1 = interp1(sigma,epsilon,x_dis);
plot(sigma,epsilon,'o','LineWidth',2);
hold on
plot(x_dis,y_interp1,'LineWidth',2);
title('interp1')

subplot(3,2,3)
y_cubicspline = cubicspline(sigma,epsilon,x_dis);
plot(sigma,epsilon,'o','LineWidth',2);
hold on;
plot(x_dis,y_cubicspline,'LineWidth',2);
title('cubicspline')

subplot(3,2,4)
y_spline = spline(sigma,epsilon,x_dis);
plot(sigma,epsilon,'o','LineWidth',2);
hold on;
plot(x_dis,y_spline,'LineWidth',2);
title('spline')

subplot(3,2,5)
y_IMQ(1,:) = polyval(polyfit(sigma,epsilon,1),x_dis);
y_IMQ(2,:) = polyval(polyfit(sigma,epsilon,2),x_dis);
y_IMQ(3,:) = polyval(polyfit(sigma,epsilon,4),x_dis);
plot(sigma,epsilon,'o','LineWidth',2);
hold on;
plot(x_dis,y_IMQ,'LineWidth',2);
title('IMQ')

subplot(3,2,6)
y_ICL = interp1(sigma,epsilon,x_dis);
plot(sigma,epsilon,'o','LineWidth',2);
hold on;
plot(x_dis,y_ICL,'LineWidth',2);
title('ICL')

%% 5
clear
clc
f = @(x) exp(-x.^2).*sin(x);
a = -2;
b = 3;
n = 3;
h = (b-a)/n;
x_dis = linspace(b,a,1000);
x_nod = a:h:b;
P_dis = interp1(x_nod,f(x_nod),x_dis);
figure()
plot(x_nod,f(x_nod),x_dis,P_dis)
err_inf = abs(f(x_dis)-P_dis);
ddf = Jac(Jac(f));
stima = @(h) ((h.^2)./8).*max(ddf(x_dis));
H = [];
err_max = [];
for n = 2.^[2:7]
    h = (b-a)/n;
    H = [H; h];
    x_nod = a:h:b;
    P_dis = interp1(x_nod,f(x_nod),x_dis);
    err_max = [err_max; max(abs(f(x_dis)-P_dis))];
end
loglog(H,err_max,'ro-',H,stima(H));
%% 6
clear
clc
y = [ 10 9.89 9.75 9.66 9.10 8.95 8.10 7.49 6.8 6.13 5.05];
f = @(t) 10 - 0.5*9.81*t.^2;
T = 0:0.1:1;
x = linspace(0,1.5,100);
n = length(y)-1;
lagr = polyval(polyfit(T,y,n),x);
int = interp1(T,y,x);
IMQ = polyval(polyfit(T,y,2),x);

subplot(2,2,1)
plot(x,f(x),'LineWidth',2)
hold on
plot(x,lagr,'LineWidth',2)
title('Lagrange')
lagr = polyval(polyfit(T,y,n),1.05)

subplot(2,2,2)
plot(x,f(x),'LineWidth',2)
hold on
plot(x,int,'LineWidth',2)
title('Interp1')
int = interp1(T,y,1.05)

subplot(2,2,3)
plot(x,f(x),'LineWidth',2)
hold on
plot(x,IMQ,'LineWidth',2)
title('IMQ')
IMQ = polyval(polyfit(T,y,2),1.05)

%% 7
clear
clc
g = @(x) 10 * x.^2;
f = @(x) g(x) + 2*rand(size(x))-1;
n = 9;
a = 0;
b = 1;
h = (b-a)/n;
x_nod = a:h:b;
x_dis = linspace(a,2,1000);
lagr = polyval(polyfit(x_nod,f(x_nod),n),x_dis);
IMQ = polyval(polyfit(x_nod,f(x_nod),2),x_dis);
plot(a:0.001:b,f(a:0.001:b), ...
    x_dis,g(x_dis), ...
    x_dis,lagr, ...
    x_dis,IMQ, LineWidth=2);
ylim([-2 20])
xlim([0 2])
lagr2 = polyval(polyfit(x_nod,f(x_nod),n),2);
IMQ = polyval(polyfit(x_nod,f(x_nod),2),2);

 