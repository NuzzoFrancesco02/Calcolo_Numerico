%% SECONDO PARZIALE
%% 1
clear
clc

x = 0:4;
y = [2 1 0 0.5 1];
spline(x,y,3.5)

%% 2
clear
clc

f = @(x) 3.*sqrt(2+abs(x));
a = -1; b = 2; h = (b-a)/3;
trapcomp(a,b,3,f)

%% 3
clear
clc
x = [-1 0 2];
phi = poly_phi(x);
expand(phi(3))

%% 4
clear
clc

n = 20;
i = 0:n;
x = 2.*(i./n)-1;
y = 3 + 4.*x + i./n;

x_dis = x(1):0.001:x(end);
y_dis = polyval(polyfit(x,y,0),x_dis);
plot(x_dis,y_dis,x,y,'o')



%% 5
clear
clc

f = @(x) x.^3.*cos(x);

H = 0.1; a = 0; b = 1; M = (b-a)/H;
trapcomp(a,b,M,f)

%% 6 
clear
clc

M1 = 10; e1 = 1e-1; e2 = 1e-6;
M2 = ceil(M1*nthroot(e1/e2,4))


%% 7
clear
clc

f = @(t,y) -2*t*y^2; y0 = 5;
tf = 10; tv = [0 tf];
[t,u] = Crank_Nicolson(f,tf,y0,0.1);
u(2)
%% 8
clear
clc

mu = 1; eta = 2; sigma = @(x) 4.*(1-x); f = @(x) 0.*x;
a = 0; b = 1; alpha = 4; beta = 0;
h = 0.1; N = 9;
[x,u] = miste(a,b,alpha,beta,mu,eta,sigma,f,N,'1','D-N','');
u(end)

%% 9
clear
clc

mu = 1; eta = 1000; sigma = @(x) 0.*x; f = @(x) 0.*x+200;
a = 0; b = 1; alpha = 5; beta = 0;
h = 1/100; N = 99;
Peh = abs(eta)*h/(2*mu);
mu = mu*(1+Peh);
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
u(100)

%% 10
clear
clc

mu = 1;
f = @(x,t) 0; a = 0; b = 1; u_s = @(t) 5; u_d = @(t) 0;
g_0 = @(x) 5.*(1-x).^2;
h = 0.5; T = 5; delta_t = 0.1;
theta = 0;
[u, x, t] = EqCalore_DiffFin_Theta(mu, f, a, b, u_s, u_d, g_0, T, h, delta_t,theta);

%% ESERCIZIO
clear
clc

f1 = @(t) -(246/25*cos(2.*t)+16.*cos(3.*t)+24/5.*sin(2.*t)).*exp(-t./5);
f2 = @(t) -(12*cos(2.*t)+1452/25*cos(3.*t)+132/5*sin(3.*t)).*exp(-t./5);
A = [-2/3 -3 0 4/3;...
    1 0 0 0;...
    0 2 -3/2 -2;...
    0 0 1 0];
g = @(t) [f1(t)./3 0 f2(t)./2 0]';
y0 = [-3/5 3 -4/5 4]'; tf = 20; tv = [0 tf];

fun = @(t,y) A*y + g(t);
h = 0.1;
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/h);
sol1 = [u(:,2),u(:,end)]
% Con un plot verifico le disuguaglianze:
plot(t,abs(u(2,:)),'-o',t,abs(u(4,:)),'-o','LineWidth',2); yline(0.5);
legend('u2','u4')
% Dal grafico si può vedere che la condizione è soffisfatta per entrambi al
% tempo tm = 19, in alternativa:
i = find(abs(u(2,:))>=0.5 | abs(u(4,:))>0.5,1,'last');
tm = t(i+1) % come già visto dal grafico

x1 = @(t) 3.*cos(2.*t).*exp(-t./5);
x2 = @(t) 4.*cos(3.*t).*exp(-t./5);
% verifico con h = 0.001
[t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/0.001);
subplot(2,1,1)
plot(t,u(2,:),t,x1(t))
subplot(2,1,2)
plot(t,u(4,:),t,x2(t))
H = [10 5 2.5 1.25].*1e-3;
Eh = [];
for i = H
    [t,u] = eulero_avanti_sistemi(fun,tv,y0,tf/i);
    Eh = [Eh max(abs(u(4,:)-x2(t)))];
end
Eh
% Se l'errore associato ad un metodo numerico per EDO è tale per cui:
%               e(n)<=C*h^p, per n = 0,...,Nh
% allora:       e(n) ~ C*h^p --> log(e(n1)/e(n2)) ~ p * log_{h(1)/h(2)}(1) 
%               --> p ~= log(e(n1)/e(n2))/log(h(1)/h(2))
%
p = log(Eh(end)/Eh(end-1))/log(H(end)/H(end-1));
% Dalla teoria, poiché f(t,y) ∈ C2([t0,tf]) in y(i) con i = 1,...,ordine del sys
% Eulero in avanti converge con ordine p = 1 in h.

% o anche:
figure()
loglog(H,Eh,H,H,'-',H,H.^2,'-.','LineWidth',2)
legend('Eh','H','H') % Si vede che Eh è parallelo ad H
%%
% In questo caso f(t,y) è differenziabile con continuità rispetto a 
% yi con i = 1,...,ord del sys e la jacobiana di f(t,y) coincide con la 
% matrice A, cerco l'autovalore massimo e il minimo di A: poiché sono entrambi
% finiti e a parte reale negativa, dalla teoria h_max è limitato 
% dall'autovalore massimo l_max. 
%
% Utilizzando la funzione di stabilità |R(z)|<1 con:
%                       z = 1 + h*l
% posso trovare h_max. In particolare per Eulero in avanti se l_max ∈ C -2*real(l_max)/abs(l_max)^2
A = [-2/3 -3 0 4/3;...
    1 0 0 0;...
    0 2 -3/2 -2;...
    0 0 1 0];
l_min = min(eig(A)); l_max = max(eig(A));

h_max = -2*real(l_max)/abs(l_max)^2
% si puo verificare anche per via grafica utilizzando R(z):
h = 0:0.001:1;
plot(h,abs(1+h.*l_max),'LineWidth',2);yline(1)
% o anche con:
options = optimoptions('fsolve','Display','off');
fsolve(@(h) abs(1+h*l_max)-1,0.4,options)
%%
clear
clc
A = [-2/3 -3 0 4/3;...
    1 0 0 0;...
    0 2 -3/2 -2;...
    0 0 1 0];
f1 = @(t) -(246/25*cos(2.*t)+16.*cos(3.*t)+24/5.*sin(2.*t)).*exp(-t./5);
f2 = @(t) -(12*cos(2.*t)+1452/25*cos(3.*t)+132/5*sin(3.*t)).*exp(-t./5);
g = @(t) [f1(t)./3 0 f2(t)./2 0]';
y0 = [-3/5 3 -4/5 4]'; tf = 20; tv = [0 tf];
h = 0.1;
x1 = @(t) 3.*cos(2.*t).*exp(-t./5);
x2 = @(t) 4.*cos(3.*t).*exp(-t./5);

[ t, u ] = eulero_indietro_sistemi_LU(A, g, tv, y0, tf/h);
sol6 = [u(4,2),u(4,3)]
% verifico con h = 0.001
[t,u] = eulero_indietro_sistemi_LU(A, g,tv,y0,tf/0.001);
figure()
subplot(2,1,1)
plot(t,u(2,:),t,x1(t))
subplot(2,1,2)
plot(t,u(4,:),t,x2(t))

[ t, u ] = multipasso_esame( A, g, tv, y0, tf/h );
sol7 = [u(4,2),u(4,3),u(4,end)]
% verifico con h = 0.001
[t,u] = multipasso_esame(A, g,tv,y0,tf/0.001);
figure()
subplot(2,1,1)
plot(t,u(2,:),t,x1(t))
subplot(2,1,2)
plot(t,u(4,:),t,x2(t))
%%
A=[-2/3 -3 0 4/3 ;1 0 0 0; 0 2 -3/2 -2; 0 0 1 0];
tf=20;
f1=@(t) -exp(-t/5).*(328/25 *cos(2*t) +20*cos(3*t) +32/5*sin(2*t));
f2=@(t) -exp(-t/5).*(16 *cos(2*t) +363/5*cos(3*t) +33*sin(3*t));
g=@(t) [1/3 * f1(t); 0; 0.5 * f2(t);0];
Afun=@(y) A+0*y;
fun=@(t,y) Afun(y)*y+g(t);
y0=[-4/5; 4; -1; 5]; tv = [0 tf]; h = 0.1;
[ t, u ] = Giulio(A, g, tv, y0, tf/h);
sol6 = [u(4,2),u(4,3)]
%%
A = [5 -1 0 0 0; 1 5 0 0 0; 0 -1 4 0 0; 0 0 1 7 0; 0 0 0 9 -1];
shift_finder(A,4)