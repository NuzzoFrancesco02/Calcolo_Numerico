%% serie08
f = @(x) exp(3.*x).*(-4+3.*x+9*x.^2);
u_ex = @(x) exp(3.*x).*(x-x.^2)+1-3.*x;
du = @(x) 3.*exp(3.*x).*(x-x.^2)-exp(3.*x).*(2.*x-1)-3;
ddu = @(x)3.*exp(3.*x).*(x-x.^2)-exp(3.*x).*(2.*x-1)-3;
a = 0;
bV = 1;
alpha = 1;
beta = -2;
mu = 1;
N = 19;
h = (bV-a)/(N+1);
x_nodes = a:h:bV;
A = sparse(1 : N, 1 : N, 2, N,N) + ...
    sparse(2 : N, 1 : N-1,-1 ,N,N) + ...
    sparse(1:N-1,2:N,-1,N,N);
A = mu/h^2*A;

bv = f(x_nodes(2:end-1))';
bv(1) = bv(1) + mu/h.^2* alpha;
bv(end) = bv(end) + mu/h.^2 * beta;
u_h = A \ bv;
u_h = [alpha; u_h; beta];
x_dis = a:0.001:bV;
plot(x_dis,u_ex(x_dis),x_nodes,u_h,'m--','LineWidth',2);
hold on;
grid on;
N = [10 20 40 80 160];
err = [];
K = [];
for i = N
    h = (bV-a)/(i+1);
    x_nodes = a:h:bV;
    A = sparse(1 : i, 1 : i, 2, i,i) + ...
        sparse(2 : i, 1 : i-1,-1 ,i,i) + ...
        sparse(1:i-1,2:i,-1,i,i);
    A = mu/h^2*A;
    
    bv = f(x_nodes(2:end-1))';
    bv(1) = bv(1) + mu/h.^2* alpha;
    bv(end) = bv(end) + mu/h.^2 * beta;
    u_h = A \ bv;
    u_h = [alpha; u_h; beta];
    err = [err; max(abs(u_ex(x_nodes)'-u_h))];
    K = [K condest(A)];
end
h = (bV-a)./(N+1);
close all
figure(1)
loglog(h,err,h,h,'--k',h,h.^2,':k','linewidth',2);
legend('Eh','N','N^2');
figure(2)
plot(N,K,'-o','LineWidth',2)

%%
clear
clc
close all
a=0; bV=1;
alpha=0; beta=1;
mu = 1e-2; eta = 1;
f = @( x ) 0 * x;
uex = @( x ) ( exp( eta / mu * x ) - 1 ) / ( exp( eta / mu ) - 1 );
       N = 20 - 1;
       h = ( bV - a ) / ( N + 1 );
       h = 0.02;
       xnodes = linspace( a, bV, N + 2 );
       A = sparse( 1 : N, 1 : N, 2, N, N ) ...
       + sparse( 2 : N, 1 : N - 1, -1, N, N ) ...
       + sparse( 1 : N - 1, 2 : N, -1, N, N );
       A = mu / h^2 * A;
       A = A + eta / ( 2 * h ) * ( sparse( 2 : N, 1 : N - 1, -1, N, N ) ...
       + sparse( 1 : N - 1, 2 : N, 1, N, N ) );
       bv = ( f( xnodes( 2 : end - 1 ) ) )';
       bv( 1 ) = bv( 1 ) + alpha * ( mu / h^2 + eta / ( 2 * h ) );
       bv( end ) = bv( end ) + beta * ( mu / h^2 - eta / ( 2 * h ) );
uh = A \ bv;
uh = [ alpha; uh; beta ];
xplot = linspace( a, bV, 1001 );
figure( 1 );
plot( xplot, uex( xplot ), '-b', xnodes, uh, ':xr' ); 
grid on; xlabel( 'x' ); ylabel( 'uex, uh ')
legend( 'u {ex}', 'uh DF centrate', 'location', 'Best' );
%% 7
clear
clc
close all
eta = 1; mu = 1;
N = 20 - 1;
h = 0.02;
Pe_h = abs(eta)*h/(2*mu);
mu_h = mu*(1 + Pe_h);
a=0; bV=1;
alpha=0; beta=1;
f = @( x ) 0 * x;
uex = @( x ) ( exp( eta / mu * x ) - 1 ) / ( exp( eta / mu ) - 1 );

       xnodes = linspace( a, bV, N + 2 );
       A = sparse( 1 : N, 1 : N, 2, N, N ) ...
       + sparse( 2 : N, 1 : N - 1, -1, N, N ) ...
       + sparse( 1 : N - 1, 2 : N, -1, N, N );
       A = mu_h / h^2 * A;
       A = A + eta / ( 2 * h ) * ( sparse( 2 : N, 1 : N - 1, -1, N, N ) ...
       + sparse( 1 : N - 1, 2 : N, 1, N, N ) );
       bv = ( f( xnodes( 2 : end - 1 ) ) )';
       bv( 1 ) = bv( 1 ) + alpha * ( mu_h / h^2 + eta / ( 2 * h ) );
       bv( end ) = bv( end ) + beta * ( mu_h / h^2 - eta / ( 2 * h ) );
uh = A \ bv;
uh = [ alpha; uh; beta ];
xplot = linspace( a, bV, 1001 );
figure( 1 );
plot( xplot, uex( xplot ), '-b', xnodes, uh, ':xr' ); 
grid on; xlabel( 'x' ); ylabel( 'uex, uh ')
legend( 'u {ex}', 'uh DF centrate', 'location', 'Best' );
%%
a = 0; bV = 1;
alpha = 1; gamma = -exp(3);
mu = 1; sigma = @(x) 0*x; eta = 0;
f = @( x ) exp( 3 * x ) .* ( - 4 + 3 * x + 9 * x.^2 );
u_ex = @( x ) exp( 3 * x ) .* ( x - x.^2 ) + 1;
N = 100 - 1;
h = ( bV - a ) / ( N + 1 );
x_nodes = linspace( a, bV, N + 2 );
A = sparse( 1 : N, 1 : N, 2, N + 1, N + 1 ) ...
+ sparse( 2 : N, 1 : N - 1, -1, N + 1, N + 1 ) ...
+ sparse( 1 : N, 2 : N + 1, -1, N + 1, N + 1 );
A = mu / h^2 * A;
A = A + eta / ( 2 * h ) * ( sparse( 2 : N, 1 : N - 1, -1, N + 1, N + 1 ) ...
+ sparse( 1 : N, 2 : N + 1, 1, N + 1, N + 1 ) );
A = A + sparse( 1 : N, 1 : N, sigma( x_nodes( 2 : end - 1 ) ), N + 1, N + 1 );
%A( end, N - 1 : N + 1 ) = mu / h * [ -1 1 ];
A( end, N - 1 : N + 1 ) = mu / (2*h) * [1 -4  3 ];
bv = f(x_nodes(2 : end))';
bv(1) = bv(1) + alpha*(mu/h^2 + eta/(2*h));
bv(end) = gamma;

u_h = A\bv;
u_h = [alpha; u_h];
x_dis = a:0.001:bV;
plot(x_dis,u_ex(x_dis),x_nodes,u_h,'--m','LineWidth',2)
%%
N = 20 - 1;
a = 0; b = 1; alpha = 1; beta = -2;
h = (b-a)/(N+1);
f = @(x) exp(3.*x).*(-4 +3.*x + 9.*x.^2);
mu = 1;
x_nodes = linspace(a,b,N+2);
A = sparse(1 : N, 1 : N, 2, N, N) + ...
    sparse(2 : N, 1 : N -1, -1, N, N) + ...
    sparse(1 : N -1, 2 : N, -1, N, N);
u_ex = @(x) exp(3.*x).*(x-x.^2)+1-3.*x;
A = A*(mu/h^2);
bV = (f(x_nodes(2:end-1))');
bV(1) = bV(1)+mu*alpha/h^2;
bV(end) = bV(end)+mu*beta/h^2;
uh = A \bV;
uh = [alpha; uh; beta];
x_dis = linspace(a,b,1000);
plot(x_dis,u_ex(x_dis),x_nodes,uh,'LineWidth',2);
legend('ESAT.','APPROS.')

NV = [10 20 40 80 160];
eh = [];
hv = [];
K = [];
for N = NV
    h = (b-a)/(N+1);
    hv = [hv h];
    x_nodes = linspace(a,b,N+2);
    A = sparse(1 : N, 1 : N, 2, N, N) + ...
        sparse(2 : N, 1 : N -1, -1, N, N) + ...
        sparse(1 : N -1, 2 : N, -1, N, N);
    u_ex = @(x) exp(3.*x).*(x-x.^2)+1-3.*x;
    A = A*(mu/h^2);
    bV = (f(x_nodes(2:end-1))');
    bV(1) = bV(1)+mu*alpha/h^2;
    bV(end) = bV(end)+mu*beta/h^2;
    uh = A \bV;
    uh = [alpha; uh; beta];
    eh = [eh max(abs(u_ex(x_nodes)-uh'))];
    K = [K condest(A)];
end
figure()
loglog(hv,eh,hv,hv,':k',hv,hv.^2,'-.k',LineWidth=2)
legend('Eh','H','H^2')
figure()
loglog(hv,K,hv,1./hv.^2,':',LineWidth=2)


N = 159;
h = (b-a)/(N+1);
x_nodes = linspace(a,b,N+2);
A = sparse(1 : N, 1 : N, 2, N, N) + ...
    sparse(2 : N, 1 : N -1, -1, N, N) + ...
    sparse(1 : N -1, 2 : N, -1, N, N);
u_ex = @(x) exp(3.*x).*(x-x.^2)+1-3.*x;
A = A*(mu/h^2);
bV = (f(x_nodes(2:end-1))');
bV(1) = bV(1)+mu*alpha/h^2;
bV(end) = bV(end)+mu*beta/h^2;
uh = A \bV;
err = max(abs(u_ex(x_nodes(2:end-1)')-uh))
condest(A);
res_norm = norm(bV-A*uh)/norm(bV);
err_rel = condest(A)*norm(uh)*res_norm
%% 1.2
clear
clc
mu = 1;
a = 0; b = 1; alpha = 1; beta = exp(3);
sigma = @(x) 1 + sin(2*pi.*x);
f = @(x) exp(3.*x).*(sin(2.*pi.*x)-8);
u_ex = @(x) exp(3.*x);
N = 20-1;
h = (b-a)/(N+1);
x_nodes = linspace(a,b,N+2);
A = sparse(1:N,1:N,2,N,N) + ...
    sparse(2:N,1:N-1,-1,N,N) + ...
    sparse(1:N-1,2:N,-1,N,N);
A = A*(mu/h^2);
A = A + sparse(diag(sigma(x_nodes(2:end-1))));
bv = (f(x_nodes(2:end-1))');
bv(1) = bv(1) + mu*alpha/h^2;
bv(end) = bv(end) + mu*beta/h^2;
uh = A\bv;
uh = [alpha;uh;beta];
x_dis = a:0.001:b;
plot(x_nodes,uh,x_dis,u_ex(x_dis),'LineWidth',2)

hv = [0.04 0.02 0.01 0.005 0.0025];
eh = [];
Nv = round((b-a)./(hv))-1;
for N = Nv
    h = (b-a)/(N+1);
    x_nodes = linspace(a,b,N+2);
    A = sparse(1:N,1:N,2,N,N) + ...
        sparse(2:N,1:N-1,-1,N,N) + ...
        sparse(1:N-1,2:N,-1,N,N);
    A = A*(mu/h^2);
    A = A + sparse(1 : N, 1 : N,sigma(x_nodes(2:end-1)),N,N);
    bv = (f(x_nodes(2:end-1))');
    bv(1) = bv(1) + mu*alpha/h^2;
    bv(end) = bv(end) + mu*beta/h^2;
    uh = A\bv;
    uh = [alpha;uh;beta];
    eh = [eh max(abs(u_ex(x_nodes)-uh'))];
end
loglog(hv,eh,hv,hv,':',hv,hv.^2,'-.','LineWidth',2)
legend('Eh','H','H^2');
%% 2.1
clear
clc

a = 0; b = 1; alpha = 0; beta = 1;
f = @(x) 0.*x;
mu = 1; eta = 1;
u_ex = @(x) (exp(eta./mu.*x)-1)./(exp(eta/mu)-1);
N = 20-1;
x_nodes = linspace(a,b,N+2);
h = (b-a)/(N+1);
A = sparse(1:N,1:N,2,N,N) + ...
    sparse(2:N,1:N-1,-1,N,N) + ...
    sparse(1:N-1,2:N,-1,N,N);
A = A*mu/h^2;
A = A + (eta/(2*h))*(sparse(2:N,1:N-1,-1,N,N) + ...
    sparse(1:N-1,2:N,1,N,N));
bv = (f(x_nodes(2:end-1)))';
bv(1) = bv(1) + alpha*(mu/h^2 + eta/(2*h));
bv(end) = bv(end) + beta*(mu/h^2 - eta/(2*h));
uh = A\bv;
uh = [alpha;uh;beta];
x_dis = a:0.001:b;
plot(x_dis,u_ex(x_dis),x_nodes,uh,'LineWidth',2);
legend('Esatta','Approx')
%% 3.1
mu = 1; eta = 0;
sigma = @(x) x.*0;
f = @(x) exp(3.*x).*(-4+3.*x+9.*x.^2);
a = 0; b = 1; alpha = 1; gamma = -exp(3);
u_ex = @( x ) exp( 3 * x ) .* ( x - x.^2 ) + 1;
N = 100-1;
h = (b-a)/(N+1);
x_nodes = linspace(a,b,N+2);
A = sparse(1:N,1:N,2,N+1,N+1) + ...
    sparse(2:N,1:N-1,-1,N+1,N+1) + ...
    sparse(1:N,2:N+1,-1,N+1,N+1);
A = mu/h^2 * A;
A = A + eta/(2*h)*(sparse(2:N,1:N-1,-1,N+1,N+1) + ...
    sparse(1:N,2:N+1,1,N+1,N+1));
A = A + sparse(1:N,1:N,sigma(x_nodes(2:end-1)),N+1,N+1);
%A(end,N:N+1) = mu/h * [-1 1];
A(end,N-1:N+1) = mu/(2*h) * [1 -4 3];
bv = f(x_nodes(2:end))';
bv(1) = bv(1) + alpha+ (mu/h^2+eta/(2*h));
bv(end) = gamma;
u_h = A\bv;
u_h = [alpha;u_h];
figure()
x_dis = a:0.001:b;
plot(x_dis,u_ex(x_dis),x_nodes,u_h,'LineWidth',2)
legend('Esatta','Approx')


%% 4
clear all
close all
clc
a = 0; b = 1; c = 0; d = 1;
uex = @(x,y) sin(pi.*x).*sin(2.*pi.*y);
h = 1e-1;
xvect = a:h:b;
yvect = c:h:d;
[Xplot,Yplot] = meshgrid(xvect,yvect);
surf(Xplot,Yplot,uex(Xplot,Yplot),'Lines','no')
xlabel('x')
ylabel('y')
zlabel('u_{ex}')
colorbar

mu = 1;
f = @(x,y) 5*pi^2*sin(pi*x).*sin(2*pi*y);
g = @(x,y) 0;
Omega = [a,b,c,d];
hx = 1e-1; hy = 1e-1;
[X,Y,U] = Poisson_Dirichlet_diff_finite_5punti(mu,f,g,Omega,hx,hy);
figure()
surf(X,Y,U)
xlabel('x')
ylabel('y')
zlabel('u_h');
%%
hvect = 0.1*2.^[0:-1:-4];
errvect = zeros(size(hvect));
for i = 1 : length(hvect)
    h = hvect(i);
    [X,Y,U] = Poisson_Dirichlet_diff_finite_5punti(mu,f,g,Omega,h,h);
    errvect(i) = max(max(abs(U-uex(X,Y))));
end
loglog(hvect,errvect,'o-','LineWidth',3)
hold on;
loglog(hvect,hvect,'k--','LineWidth',3)
loglog(hvect,hvect.^2,'k-.','LineWidth',3)
grid on
legend('errore','lineare','quadratico')
%% 
a = 0; b = 10; 
T = 5;
mu = 1;
f = @(x,t) (-sin(t) + 0.25*cos(t))*sin(0.5*x);
u_ex = @(x,t) sin(0.5*x)*cos(t);
u_s = @(t) 0*cos(t);
u_d = @(t) sin(0.5)*cos(t);
g_0 = @(x) sin(0.5*x);
theta = 1; h = 0.01; delta_t = 0.01;
[u,x,t] = Merda(mu, f, a, b, u_s, u_d, g_0, T, h, delta_t,theta);
figure;
for i = 2 :length(t)
    plot(x,u(:,i),'-',x,u_ex(x,t(i)),'LineWidth',3)
    pause(t(i)-t(i-1));
    axis([0 10 -1 1])
end
legend('APPROX','ESATT')
title(['Eulero Implicito, h = ',num2str(h),' \Delta_t = ',num2str(delta_t)]);
%%
a = 0; b = 1; alpha = 0; beta = 0; gamma = 0; b_b = 3;
eta = 0; sigma = @(x) 0.*x; f = @(x) 6 + 0.*x; N = 9; mu=1;
[x_nodes,uh] = Robin(a,b,alpha,gamma,b_b,mu,eta,sigma,f,N,1e-1) ;
plot(x_nodes,uh)
%%
sigma = @(x) 0.*x;
f = @(x) 4 + 0.*x;
[x_nodes,u] = Robin(0,1,0,2,3,1,0,sigma,f,9,1/10);
u(2)
