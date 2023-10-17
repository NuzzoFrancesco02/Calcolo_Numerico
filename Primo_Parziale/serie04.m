%% serie04
%% 1
f = @(x) x.^3 - (2+exp(1)).*x.^2 +(2*exp(1)+1).*x+(1-exp(1))-cosh(x-1);
x = 0.5:0.01:6.5;
plot(x,f(x));
grid on;
zero_1 = bisez(3,4,1e-12,f);
zero_2 = bisez(6,6.5,1e-12,f);

%% 2
f = @(x) x.^3 - (2+exp(1)).*x.^2 +(2*exp(1)+1).*x+(1-exp(1))-cosh(x-1);
df = @(x) 3*x.^2 - 2*(2+exp(1)).*x + (2*exp(1)+1) - sinh(x-1);
x = 0.5:0.01:6.5;
plot(x,f(x),'r');
hold on;
plot(x,df(x),'b');
grid on;

[sol_1 , iter]= newton(0.5,1e-6,1000,f,df);
[sol_1_m, iter_m] = newton(0.5,1e-6,1000,f,df,2);
sol_2 = newton(3,1e-6,1000,f,df);
sol_3 = newton(6,1e-6,1000,f,df);

alpha1 = 1;
err_m = abs(sol_1_m-alpha1);
err = abs(sol_1-alpha1);
figure(2)
semilogy(0 : iter, err, 'b', 0 : iter_m, err_m, 'r');
grid on;

%% 3
% syms f(x)
% f(x) = atan(7.*(x-pi/2)) + sin((x-pi/2).^3);
% simplify(diff(f(x)))

f = @(x) atan(7.*(x-pi/2)) + sin((x-pi/2).^3);
df = @(x) 7./((7.*x - (7.*pi)./2).^2 + 1) + 3.*cos((x - pi/2).^3).*(x - pi/2).^2;
x = -1:0.01:6;
subplot(2,2,1)
plot(x,f(x),'r',x,df(x),'b');
grid on;
[sol_1, iter_1] = newton(1.5,1e-10,1000,f,df);
[sol_2, iter_2] = newton(4,1e-10,1000,f,df);
err1 = abs(sol_1-(pi/2));
err2 = abs(sol_2-(pi/2));
subplot(2,2,2)
plot(0:iter_1,err1,'r');
subplot(2,2,3)
plot(0:iter_2,err2,'b');
a = -1;
b = 6;
[sol_bi, iter_bi]= bisez(b,a,(b-a)/(2^31),f);
err_bi = abs(sol_bi-(pi/2));
subplot(2,2,4)
plot(0:iter_bi,err_bi,'g');
[sol_comb, iter_comb] = biseznewton(a,b,5,1000,1e-10,f,df);
err_combi = abs(sol_comb-(pi/2));
figure(2)
plot(0:iter_comb,err_combi);
%% 4
f = @(x) cos(2.*x).^2 - x.^2;
df = @(x) -4.*cos(2.*x).*sin(2.*x) - 2.*x;
x = -4:0.01:4;
figure(1)
plot(x,f(x))
grid on;
% |phi(alpha)'|<1   |1 + A*f(alpha)'| < 1   (-2/f(alpha)' < A < 0)
% poiché pendenza negativa nell'alpha considerato, f(alpha)' è negativo: si
% ribaltano le disequazioni:
% 0 < A < -2/f(alpha)'
A = 0.1;
phi = @(x) x + A .* f(x);
figure(2)
subplot(3,2,1);

[succ1,it1] = ptofis(0.1,phi,1000,1e-10,-pi/2,pi/2);
title(['A = 0.1 x0 = 0.1, iter = ',num2str(it1)]);
hold on;
sup = -2/df(succ1(end));
inf = 0;
fprintf('%f < A < %f \n', inf, sup);
subplot(3,2,2);

phi = @(x) x + 0.6 .* f(x);
[succ2,it1] = ptofis(0.1,phi,1000,1e-10,-pi/2,pi/2);
title(['A = 0.6 x0 = 0.1, iter = ',num2str(it1)]);
subplot(3,2,3);
phi = @(x) x + 0.6 .* f(x);
[succ3,it1] = ptofis(2,phi,1000,1e-10,-pi/2,pi/2);
title(['A = 0.6 x0 = 2, iter = ',num2str(it1)]);
subplot(3,2,4);
phi = @(x) x + 0.8 .* f(x);
[succ4,it1] = ptofis(0.1,phi,1000,1e-10,-pi/2,pi/2);
title(['A = 0.8 x0 = 0.1, iter = ',num2str(it1)]);
subplot(3,2,5)
phi = @(x) x + 0.715426 .* f(x);
[succ5,it1] = ptofis(0.1,phi,1000,1e-10,-pi/2,pi/2);
title(['A = 0.715426 x0 = 0.1, iter = ',num2str(it1)]);

[p1 c1] = stimap(succ1);
[p2 c2] = stimap(succ2);
[p3 c3] = stimap(succ3);
[p4 c4] = stimap(succ4);
[p5 c5] = stimap(succ5);
p1(end)
c1(end)

p2(end)
c2(end)

p3(end)
c3(end)

p4(end)
c4(end)

p5(end)
c5(end)
%%
f = @(x) cos(2.*x).^2 - x.^2;
df = @(x) -4.*cos(2.*x).*sin(2.*x) - 2.*x;

A_opt = -1/df(succ1(end));
phi_opt = @(x) x + A_opt .* f(x);
succ1 = ptofis(0.1,phi_opt,1000,1e-10,-pi/2,pi/2);
stimap(succ1);
%%
phi_newton = @(x) x - f(x)/df(x);

