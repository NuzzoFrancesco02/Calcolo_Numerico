%% AUTOVALUTAZIONE 8
%% 1
f = @(x) 1 - exp(sin(x));
x_nod = linspace(-2,2,1000);
[y_pol,~, err] = poly_lagr(f,5,-2,2,x_nod,'equi');
[err_max, i] = max(err);

x_nod(i)
err_max

%% 2
n = 3;
a = 0;
b = 1;
str = 'CGL';
[~,~,P] = poly_equation(n,a,b,str) 
%% 3
f = @(x) cos(pi.*x);
a = 0;
b = 1;
n = 3
k = 0:n;
t = -cos(pi*k/n);
x = ((b-a)/2)*t + (a+b)/2;
[~,phi] = poly_equation(x)
P = simplify(expand(phi));
phi0 = matlabFunction(P(1));
phi1 = matlabFunction(P(2));
phi2 = matlabFunction(P(3));
phi3 = matlabFunction(P(4));
x = a:0.01:b;
max(abs(phi0(x))+abs(phi1(x))+abs(phi2(x))+abs(phi3(x)))
%% 2
f = @(x) exp(3.*x);
poly_stima(f,3,-1,1,'CGL')
x = -1:0.001:1;
y = lagr_polin(3,f,-1,1,x,'CGL');
err = max(abs(f(x)-y))
%% 3
A = lebesgue_CGL(3,1)
%% 4
x_nod = [997 1975 2153 2722 2431 2546 1782 1040 1670];
media_mobile_semplice(x_nod,7,5)
%% 5
f = @(x) exp(x).*sin(pi.*x);
x_nod = linspace(-1,1,11);
x = linspace(-1,1,1000);
y = interp1(x_nod,f(x_nod),0.7);
[y,~,y0]=interp_tratti(f,-1,1,1,10,'equi',0.7);


%% 6
f = @(x) (x./pi).*(x<=pi) + ((2-x/pi).^2).*(x>pi);
I_f = poly_trigo(f,4)
%x = linspace(0,2*pi,1000);
%plot(x,f(x),x,I_f(x))
%% 7
cost = 1e-2; 
M1 = 10;
M2 = ceil(nthroot((cost*M1^4/1e-5),4));
e1 = 1e-2;
e2 = 1e-5;
M2 = ceil(M1*nthroot(e1/e2,4))
%%
x = 0 : 0.001: 4*pi;
plot(x,cos(pi.*x))
%% 1
f = @(x) 1 -exp(sin(x));
a = -2;
b = 2;
n = 5;
h = (b-a)/n;
x_nod = a:h:b;
P = polyfit(x_nod,f(x_nod),n);
x = -2:0.001:2;
y_dis = polyval(P,x);
err = abs(f(x)-y_dis);
[err_max, i] = max(err);
err_max
x(i)

%% 6
f = @(x) (x./pi).*(x<=pi) + ((2-x./pi).^2).*(x>pi);
trigonometric(f,4)
%% 7 
%% 1
a = -2;
b = 2;
f = @(x) 1-exp(sin(x));
x = a:0.00001:b;
y = poly_lagr(f,5,a,b,x,'equi');
[err_max,i] = max(abs(f(x)-y));
err_max
x(i)
%% 3
[~,~,A] = poly_equation(3,0,1,'CGL');
A = sym(A)
%% 4
x = [997 1975 2153 2722 2431 2546 1782 1040 1670];
T = 7;
t = 5;
media_mobile_semplice(x,T,t);
T = 7;
t = 5;
n = (T-1)/2;
S = 0;

for j = -n : n
    S = S + x(t+j);
end
x_tilde = S/(2*n+1)
%% 5
f = @(x) exp(x).*sin(pi.*x);
a = -1;
b = 1;
h = (b-a)/11;
x_nod = a:h:b;
y0 = interp1(x_nod,f(x_nod),0.7)
%% 6
f = @(x) (x./pi).*(x<=pi) + ((2-x./pi).^2).*(x>pi);
poly_trigo(f,4)

%% 7
M1 = 10;
ceil(M1*nthroot(1e-2/1e-5,4))
%% 8
x = linspace(0,1,150);
rng(1);
y = 3*x.^2+0.3*sin(100*pi*x)+0.3*randn(1,150);
polyval(polyfit(x,y,2),1.5)
%%
n = 3; i = 0:3; xnodes = 0.5 + 0.5*(-cos(pi*i/3));
a1 = polyfit(xnodes,[1 0 0 0],3);
a2 = polyfit(xnodes,[0 1 0 0],3);
a3 = polyfit(xnodes,[0 0 1 0],3);
a4 = polyfit(xnodes,[0 0 0 1],3);
xval = linspace(0,1,1000);
max(abs(polyval(a1,xval)) + abs(polyval(a2,xval))+abs(polyval(a3,xval))+abs(polyval(a4,xval)))

%%
n = 3; a = 0; b = 1; str = 'CGL';
[eq, phi, A] = poly_equation(n,a,b,str);
A
%%
M1 = 10; e1 = 1e-2; e2 = 1e-5;
M2 = M1*nthroot(e1/e2,4)
%% 
hl = -4 : 0.001 : 4;
R = @(hl) 1 + hl + hl.^2./2 + hl.^3./6 + hl.^4./24;
plot(hl,abs(R(hl))); yline(1)
