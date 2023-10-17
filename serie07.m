%% serie07
%% 1 
f = @(x) exp(-x.^2).*sin(2*x+1);
df = Jac(f);

xb = 0;
dfxb = df(xb);
h_v = 0.4*2.^(-[0:5]);
dfa_v = [];
dfi_v = [];
dfc_v = [];
for h = h_v
    dfa_v = [dfa_v, (f(xb+h)-f(xb))/h];
    dfi_v = [dfi_v, (f(xb)-f(xb-h))/h];
    dfc_v = [dfc_v, (f(xb+h)-f(xb-h))/(2*h)];
end
err_dfa_v = abs(dfa_v-dfxb);
err_dfi_v = abs(dfi_v-dfxb);
err_dfc_v = abs(dfc_v-dfxb);

loglog(h_v,err_dfa_v,'-ob',h_v,err_dfi_v,'-xr',...
    h_v,err_dfc_v,'-sg',h_v,h_v,'--k',h_v,(h_v).^2,'-.k')
%% 2
lambda = 2.4;
f = @(t,y) lambda*y;
y0 = 0.1;
t_max = 3;

f = @(t) y0*exp(lambda*t);
figure()
t_plot = 0:0.01:t_max;
plot(t_plot,f(t_plot),'k')

h = 0.05;
[t_hEI1,u_hEI1] = eulero_avanti(f,t_max,y0,h)
h = 0.01;
[t_hEI2,u_hEA2] = eulero_avanti(f,t_max,y0,h)

hold on
plot(t_hEI1,u_hEI1,'o');
plot(t_hEI2,u_hEA2,'or');
xlabel('t');
ylabel('y');
legend('y(t)','EA h = 0.05','EA h = 0.01')

%%
lambda = 2.4;
f = @(t,y) lambda*y;
y0 = 0.1;
t_max = 3;

f = @(t) y0*exp(lambda*t);
figure()
t_plot = 0:0.01:t_max;
plot(t_plot,f(t_plot),'k')

h = 0.05;
[t_hEI1,u_hEI1,iter_pt1] = eulero_indietro(f,t_max,y0,h);
h = 0.01;
[t_hEI2,u_hEI2,iter_pt2] = eulero_indietro(f,t_max,y0,h);

hold on
plot(t_hEI1,u_hEI1,'o');
plot(t_hEI2,u_hEI2,'or');
xlabel('t');
ylabel('y');
legend('y(t)','EI h = 0.05','EI h = 0.01')

%%
figure()
plot(t_hEI1,iter_pt1,'o:')
N = 5;
dimezz = 2.^[1:N];
passi = 0.08./dimezz;

for it = 1:N
    [EA_t_h,EA_u_h] = eulero_avanti(f,t_max,y0,passi(it));
    [EI_t_h,EI_u_h,iter_pt] = eulero_indietro(f,t_max,y0,passi(it));
    y_h = f(0:passi(it):t_max);
    errore_EA(it) = max(abs(y_h-EA_u_h));
    errore_EI(it) = max(abs(y_h-EI_u_h));
end
figure()
loglog(passi,errore_EA,'-ob')
hold on;
loglog(passi,errore_EI,'-or')
loglog(passi,100*passi,'k',passi,(10*passi).^2,'k:')
legend('errore_EA','errore_EI','ordine 1','ordine 2')
%%
lambda = -42;
t_max = 1;
t_plot = 0:0.001:t_max;
y0 = 2;
f = @(t) y0*exp(lambda*t);

figure()
plot(t_plot,f(t_plot),'k','LineWidth',2)

h = 0.01;
f = @(t,y) lambda*y;
[t_hEA,u_hEA] = eulero_avanti(f,t_max,y0,h);
hold on
plot(t_hEA,u_hEA,'-or','LineWidth',2)

df = @(t,y) lambda;
[t_hEI,u_hEI,iter_new] = eulero_indietro_newton(f,df,t_max,y0,h);
plot(t_hEI,u_hEI,'-ob','LineWidth',2)

%%
y0 = 2;
lambda = -42;
t_max = 1;
t_plot = 0:0.001:t_max;
f = @(t) y0*exp(lambda*t);
f = @(t,y) lambda*y;
h = 0.02;

[t_hCN,u_hCN,iter_CN] = Crank_Nicolson(f,t_max,y0,h);

figure()
plot(t_plot,f(t_plot),'k',t_hCN,u_hCN,'ro','LineWidth',2);
legend('y(t)','u_CN')
[t_hH,u_hH] = Heun(f,t_max,y0,h);
figure()
plot(t_plot,f(t_plot),'k',t_hH,u_hH,'ro','LineWidth',2);
legend('y(t)','u_H')
%%
N = 5;
dimezz = 2.^(1:N);
H = 0.04./dimezz;
err_EI = [];
err_CN = [];
err_H = [];
err_RK = [];
for h = H
    [t_hEI,u_hEI] = eulero_indietro_pto_fisso(f,t_max,y0,h);
    [t_hCN,u_hCN,iter_CN] = Crank_Nicolson(f,t_max,y0,h);
    [t_hH,u_hH] = Heun(f,t_max,y0,h);
    [t_hRK,u_hRK] = Runge_Kutta_4(f,t_max,y0,h);
    err_EI = [err_EI, max(abs(f(t_hEI)-u_hEI))];
    err_CN = [err_CN, max(abs(f(t_hCN)-u_hCN))];
    err_H = [err_H, max(abs(f(t_hH)-u_hH))];
    err_RK = [err_RK, max(abs(f(t_hRK)-u_hRK))];
end
figure()
loglog(H,err_EI,'-o',H,err_CN,'-o',H,err_H,'-o',H,err_RK,'-o',H,H,'--',H,H.^2,'-.',H,H.^4,'LineWidth',2)
legend('EI','CN','Heun','RK','H','H^2','H^4')

%% 1
f = @(x) exp(-x.^2).*sin(2.*x+1);
h = 0.4.*(2.^(-(0:5)));
x = 0;
df = Jac(f);
ddf = Jac(df);
j = 1;
for i = h
    df_av(j) = (f(x+i)-f(x))/i;
    err_av(j) = abs(df(x)-df_av(j));
    df_ind(j) = (f(x)-f(x-i))/i;
    err_ind(j) = abs(df(x)-df_ind(j));
    df_cen(j) = (f(x+i)-f(x-i))/(2*i);
    err_cen(j) = abs(df(x)-df_cen(j));
    ddf_cen(j) = (f(x+i)-2*f(x)+f(x-i))/(i^2);
    err2_cen(j) = abs(ddf_cen(j)-ddf(x));
    j = j + 1;
end
loglog(h,err_av,'r',h,err_ind,'b',h,err_cen,'y',h,err2_cen,'g',h,h,':k',h,h.^2,'-.k')
legend('Avanti','Indietro','Centrate','Der^2 Centrate','H','H^2');
[p,c] = stimap_2(err2_cen,h);
%% 2
lambda = 2.4;
f = @(t,y) lambda*y;
t0 = 0;
t_max = 3;
y0 = 0.1;
h1 = 0.05;
h2 = 0.01;
[t1,u1]= eulero_avanti(f,t_max,y0,h1);
[t2,u2] = eulero_avanti(f,t_max,y0,h2);
es = @(t) y0.*exp(lambda*(t-t0));
t = t0:0.01:t_max;
plot(t1,u1,t2,u2,t,es(t),'LineWidth',2);
legend('EA h = 0.05','EA h = 0.01','Sol esatta');
[t3,u3,it1]= eulero_indietro(f,t_max,y0,h1);
[t4,u4,it2]= eulero_indietro(f,t_max,y0,h2);
figure()
plot(t3,u3,t4,u4,t,es(t),'LineWidth',2)
legend('EI h = 0.05','EI h = 0.01','Sol esatta');
figure()
plot(t3,it1,'o',t4,it2,'*','LineWidth',2);
legend('It h = 0.05','It h = 0.01')
%% 1.6
lambda = 2.4;
f = @(t,y) lambda*y;
t0 = 0;
t_max = 3;
y0 = 0.1;
h = 0.04./(2.^(0:4));
err_EA = [];
err_EI = [];
es = @(t) y0.*exp(lambda*(t-t0));
for i = h
    [t1,u_EA] = eulero_avanti(f,t_max,y0,i);
    err_EA = [err_EA max(abs(es(t1)-u_EA))];
    [t2,u_EI] = eulero_indietro(f,t_max,y0,i);
    err_EI = [err_EI max(abs(es(t2)-u_EI))];
end
loglog(h,err_EA,h,err_EI,'LineWidth',2)
legend('Errore EA','Errore EI')
%% 2.3
lambda = -42;
y0 = 2;
t0 = 0;
es = @(t) y0*exp(lambda*t);
f = @(t,y) lambda*y;
df = @(t,y) lambda;
t_max = 1;
h = 0.05;
t = t0:0.01:t_max;

[t_EA,u_EA] = eulero_avanti(f,t_max,y0,h);
[t_EI,u_EI] = eulero_indietro_newton(f,df,t_max,y0,h);
plot(t_EA,u_EA,'o',t_EI,u_EI,'+',t,es(t),'LineWidth',2)
legend('EA','EI')
h = 0.01;
[t_EA,u_EA] = eulero_avanti(f,t_max,y0,h);
[t_EI,u_EI] = eulero_indietro_newton(f,df,t_max,y0,h);
figure()
plot(t_EA,u_EA,'o',t_EI,u_EI,'+',t,es(t),'LineWidth',2)
legend('EA','EI')
%% 2.3.b
lambda = -42;
y0 = 2;
t0 = 0;
es = @(t) y0*exp(lambda*t);
f = @(t,y) lambda*y;
df = @(t,y) lambda;
t_max = 1;
h = 0.05;
t = t0:0.01:t_max;
[t_EA,u_EA] = eulero_avanti(f,t_max,y0,h);
[t_EI,u_EI] = eulero_indietro_pto_fisso(f,t_max,y0,h);
figure()
plot(t_EA,u_EA,'-o',t_EI,u_EI,'-+',t,es(t),'LineWidth',2);
legend('EA h = 0.05', 'EI h = 0.05')
title(num2str(lambda*h))
ylim([-2 10])

h = 0.03;
[t_EA,u_EA] = eulero_avanti(f,t_max,y0,h);
[t_EI,u_EI] = eulero_indietro_pto_fisso(f,t_max,y0,h);
figure()
plot(t_EA,u_EA,'-o',t_EI,u_EI,'-+',t,es(t),'LineWidth',2);
legend('EA h = 0.03', 'EI h = 0.03')
title(num2str(lambda*h))
ylim([-2 10])

h = 0.01;
[t_EA,u_EA] = eulero_avanti(f,t_max,y0,h);
[t_EI,u_EI] = eulero_indietro_pto_fisso(f,t_max,y0,h);
figure()
plot(t_EA,u_EA,'-o',t_EI,u_EI,'-+',t,es(t),'LineWidth',2);
legend('EA h = 0.01', 'EI h = 0.01')
title(num2str(lambda*h))
%%
lambda = -42;
y0 = 2;
t0 = 0;
es = @(t) y0*exp(lambda*t);
f = @(t,y) lambda*y;
df = @(t,y) lambda;
t_max = 1;
h = 0.02;
t = t0:0.01:t_max;
dimezz = 2.^(0:4);
h = 0.02./dimezz;
err_CN = [];
err_EI = [];
err_He = [];
err_RK = [];
for i = h
    [t_EI,u_EI] = eulero_indietro_pto_fisso(f,t_max,y0,i);
    err_EI = [err_EI max(abs(es(t_EI)-u_EI))];
    [t_CN,u_CN,it_CN] = Crank_Nicolson(f,t_max,y0,i);
    err_CN = [err_CN max(abs(es(t_CN)-u_CN))];
    [t_He,u_He] = Heun(f,t_max,y0,i);
    err_He = [err_He max(abs(es(t_He)-u_He))];
    [t_RK,u_RK] = Runge_Kutta_4(f,t_max,y0,i);
    err_RK = [err_RK max(abs(es(t_RK)-u_RK))];
end
loglog(h,err_EI,h,err_CN,h,err_He,h,err_RK,h,h,':k',h,h.^2,'-k',h,h.^4,'-.k','LineWidth',2);
legend('EI','Clark-Nicolson','Heun','Runge-Kutta','H','H^2','H^4')
%%
y0 = [2 0];
[t,u] = ode45('fmm',[0 100],y0);
% u = [x,x']
plot(t,u(:,1))
%%
t0 = 0; tf = 5;
y0 = [1 0 0 -5.1];
options = odeset('reltol',1e-6);
[tv,uv] = ode45('twobody',[t0 tf],y0,options);
uv = uv';
figure(1)

plot(tv,uv(1,:),'--k',tv,uv(3,:),'--m','LineWidth',2,'MarkerSize',10);
axis equal;
grid on;
figure()
for i = 1 : length(tv)
plot(0,0,'xk',uv(1,i),uv(3,i),'og', ...
    uv(1,:),uv(3,:),'-b','LineWidth',2,'MarkerSize',10);
pause(tv(i+1)-tv(i));
end
%figure()
%plot(tv(1:end-1),tv(2:end)-tv(1:end-1),'LineWidth',2);

