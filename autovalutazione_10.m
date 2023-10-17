%% Autovalutazione_10
%% 1
f = @(x) 5+2.^x;
x0 = 0;
h = 1/4;
df = @(x) (-3*f(x)+4*f(x+h)-f(x+2*h))/(2*h);
df(x0)
%% 2
syms g;
f = @(x) g.*(x.^3) + 7.*x - 53;
h = 0.5;
d3 = 6*g;
er  = -(h^2)*(2*d3)/12;
%% 3
syms h t y;
f = -y+2*t;
u3 = EA_sym(f,3,1)
%% 4
syms h t y;
f = -y + 2*t;
u3 = EAI_sym(f,3,1)
%%
f = @(t,y) -(1+t.*0.5).*y;
h = 0.1;
[~,u]= eulero_indietro_pto_fisso(f,0.1,3,0.1)
%%
f = @(t,y) -y.*exp(t.*y);
y0 = 2;
n = 3;
h = 0.1;
t = h.*[0:n];
u(1) = y0;
for i = 1:n
    u_mez = u(i)+0.5*h*f(t(i),u(i));
    q = f(t(i)+0.5*h,u_mez);
    u(i+1) = u(i)+h*q;
end
u(end)