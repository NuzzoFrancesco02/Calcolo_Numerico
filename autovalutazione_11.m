%% autovalutazione_11
%% 2
A = [0 0; 3/4 0];
b = [1/3 2/3];
c = [0 3/4];
f = @(t,y) -2*(1+sin(pi.*t)).*((2-y).^2);
t_max = 5;
y0 = 3;
h = 1/4;
u = Runge_Kutta_gen(A,b,c,f,t_max,h,y0);
u(end)
%% 3

A = [0 0 0 0; 1/3 0 0 0; -1/3 1 0 0; 1 -1 1 0];
b = [1/8 3/8 3/8 1/8];
c = [0 1/3 2/3 1];
Runge_Kutta_gen(A,b,c)

%%

t0 = 0;
tf = 3;
y0 = 2;
f = @(t,y) sin(pi.*t)-y.^2;
h = 0.2;

[t,u] = Adams_Bashforth(f,[0 3],h,y0,1.2,0.9120);
plot(t,u);
%%
clear
clc

t0 = 0;
tf = 10;
y0 = 1;
w0 = 0;
h = 0.001;
f = @(t,y,v) sin(pi.*t)-y.^2-3*v;
[t,u,v] = Leap_Frog(f,[0 10],0.2,1,0);
u(end)
v(end)