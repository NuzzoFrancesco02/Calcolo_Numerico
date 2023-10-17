%% AUTOVALUTAZIONE 12
%% 1
A = [-4 1; 2 -3];
y0 = [7;2];
fun = @(t,y) A*y;
tf = 1; tv = [0 tf];
[t,u] = eulero_indietro_sistemi(fun,tv,y0,tf/0.1);
u(:,end)
%% 2
f = @(t,y) [-exp(-t)*y(1)+y(2) -y(2)+t y(1)-2*y(3)]';
tf = 1; tv = [0 tf];
y0 = ones(3,1);
[t,u] = Crank_Nicolson_sistemi(f,tv,y0,tf/0.1);
u(:,end)
%% 3
A = [-2 0 0; 0 -3 1; 0 -1 -3];
l = max(eig(A));
R_abs = @(h) (abs(1 + h*l)-1);
fsolve(R_abs,1)
%% 4
mu = 2; f = @(x) 10.*x.^2; sigma = @(x) 0.*x;
eta = 0; a = 0; b = 1; alpha = 0; beta = 4;
[x,uh] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,9,0.1);
uh(10)
%% 5
mu = 1; eta = 0; sigma = @(x) 0.*x + pi^2; f = @(x)(x+2*sin(pi*x))*pi^2;
a = 0; b = 1; alpha = 1 + pi; gamma = 1 - pi;
[x,u] = miste(a,b,alpha,gamma,mu,eta,sigma,f,9,'1','N-N','');
u(end)
plot(x,u)
%% 6
mu = 3; eta = -100;
2*mu/abs(eta)
%% 7
mu = 1; eta = -100;
a = -1; b = 1; alpha = 7; beta = 0;
sigma = @(x) 0.*x; f = @(x) 0.*x;
h = 0.2;
N = (b-a)/h-1;
Peh = abs(eta)*h/(2*mu);
mu = mu*(1+Peh);
[x,u] = Dirichlet(a,b,alpha,beta,mu,eta,sigma,f,N,h);
u(2)
%% 8 
e1 = 4*1e-2; h1 = 0.1; h2 = 0.05;
e2 = e1*(h2/h1)^2
%% 9
mu = 1; f = @(t,x) 0;
u_s = @(t) 0; u_d = @(t) 0;
g_0 = @(x) 7*sin(pi*x);
T = 1; h = 0.5; delta_t = 0.2;
a = 0; b = 1;
[u, x, t] = EqCalore_DiffFin_Theta(mu, f, a, b, u_s, u_d, g_0, T, h, delta_t,1);
u(2,end)
