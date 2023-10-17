clc
clear 
close all

%% Esercizio 1 - Ricerca degli zeri con il Metodo di Newton

% Punto 1

f = @(x) x.^3 - (2+exp(1))*x.^2 + (2*exp(1)+1)*x + (1-exp(1)) - cosh(x-1);
x = linspace(0.5, 6.5, 100);
y = f(x);
figure(1);
plot(x,y)
title('f(x)=x^3 - (2+e)x^2 + (2e+1)x + (1-e) - cosh(x-1)');
xlabel('x');
ylabel('y');
grid on
y0 = zeros(100,1);
hold on
plot(x,y0)

df = @(x) 3*x.^2 - 2*(2+exp(1))*x + (2*exp(1)+1) - sinh(x-1);
dy = df(x);
d2f = @(x) 6*x - 2*(2+exp(1)) - cosh(x-1);
d2y = d2f(x);
figure(2)
plot(x, d2y, 'k', x,dy,'b', x,y, 'r', x,y0, 'g')
title('df(x) e d2f(x)');
xlabel('x');
ylabel('y');
legend('y=d2f(x)', 'y=df(x)', 'y=f(x)', 'y=0', 'Location', 'SouthWest')
grid on

% discussione sulle radici dal grafico di f e df
alpha1 = 1;
if (df(alpha1) == 0)
    if (d2f(alpha1) == 0)
        disp('la radice alpha1 ha molteplicita'' maggiore di due')
    else
        disp('la radice alpha1 ha molteplicita uguale a due')
    end
else
    disp('la radice alpha1 ha molteplicita'' uguale ad uno')
end

%% Punto 2
% vedi newtom.m

%% Punto 3

toll = 1e-6;
nmax = 100;

disp('Calcolo della radice con molteplicita'' 2')
x01 = 0.5;
disp('Metodo di Newton semplice');
[xvect_1,it_1]=newton(x01,nmax,toll,f,df);
disp('Metodo di Newton modificato');
[xvect_1m,it_1m]=newton(x01,nmax,toll,f,df,2);


disp('Calcolo della prima radice con molteplicita'' 1')
x02 = 3;
[xvect_2,it_2]=newton(x02,nmax,toll,f,df);

disp('Calcolo della prima radice con molteplicita'' 1')
x03 = 6;
[xvect_3,it_3]=newton(x03,nmax,toll,f,df);

alpha1 = 1; % soluzione esatta 

% valutazione errore con Newton e Newton modificato
sol_1 = alpha1*ones(length(xvect_1),1);
sol_1m = alpha1*ones(length(xvect_1m),1);

err_1 = abs(xvect_1 - sol_1);
err_1m = abs(xvect_1m - sol_1m);

it_1vec = 0:it_1;
it_1mvec = 0:it_1m;

figure(3);
semilogy(it_1vec,err_1,'b', it_1mvec,err_1m,'r', 'LineWidth', 2);
grid on   
title('Confronto convergenza metodo Newton e Newton modificato')
legend('Newton','Newton modificato','Location', 'NorthEast');
xlabel('iterazioni');
ylabel('errore');