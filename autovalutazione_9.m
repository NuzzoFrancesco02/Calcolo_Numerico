%% AUTOVALUTAZIONE 9
%% 1
x = linspace( 0, 1, 150 );
rng( 1 );
y = 3 * x.^2 + 0.3 * sin( 100 * pi * x ) + 0.3 * randn( 1, 150 );
y0 = polyval(polyfit(x,y,2),1.5)
%% 2
f = @(x) cos(pi.*x);
I = simpcomp(-1,1,1,f)
%% 3
f = @(x) exp(-x);
M = ceil(sqrt(((2^3)/24)*f(-1)/(1e-2)))
%% 4
f = @(x) exp(x) + exp(sin(2.*pi.*x));
a = -2;
b = -2;
I = gausslegendre_comp(-2,2,f,8,1,'equi')
%% 5
f = @(x) 1+x;
clenshaw_curtis(f,2)
%% 1
x = linspace(0,1,150);
rng(1);
y = 3*x.^2 + 0.3*sin(100*pi*x)+0.3*randn(1,150);
polyval(polyfit(x,y,2),1.5);
%% 2
f = @(x) cos(pi.*x);
I = simpcomp(-1,1,1,f)
%% 3
syms x beta g
f =  exp(-x)-beta*x+g;
b = 1;
a = -1;
integr_stima(-1,1,1e-2,f,'pm')
%% 4
f = @(x) exp(x) + exp(sin(2*pi.*x));
gausscomp(-2,2,8,f)
%% 5
n = 2;
f = @(x) 1+x;
for k = 0 : n
    fun = @(t) f(cos(t)).*cos(k*t);
    a(k+1) = 2*simpcomp(0,pi,1e5,fun)/pi;
end
S = a(1)
for k = 1 : (n/2)
    S = S + a(2*(k)+1)/(1-(2*k))
end
S
%%
L = tf(0.5.*[1 1],conv([0.1 1],[0.1 1]));
bode(L)
figure()
nyquist(L)
%% 1
f = @(x) cos(pi.*x);
I = simpcomp(-1,1,1,f);
I = sym(I)
%% 2
b = 1;
a = -1;
f = @(x) exp(-x)-x+1;
ddf = @(x) exp(-x);
toll = 1e-2;
N = ceil(sqrt(((b-a)^3)*max(abs(ddf(x)))/(24*toll)))
%% 3
f = @(x) exp(x) + exp(sin(2.*x.*pi));
gausslegendre_comp(-2,2,f,8,1,'equi')
%% 4
n = 2;
k = 0:n;
f = @(x) (1+x);
a = [];
for k = 0:n
    fun = @(t) f(cos(t)).*cos(k.*t);
    a = [a 2*gausslegendre_comp(0,pi,fun,1000,9,'equi')/pi];
end
S = a(1);
for k = 2:n/2-1
    S = S + a(2*k)/(1-(2*k)^2);
end
S
clenshaw_curtis(f,2)
