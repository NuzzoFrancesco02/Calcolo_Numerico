clc
clear
close all

% Esercizio 1 - Ricerca degli zeri con il Metodo di Bisezione

%% Punto 1

f = @(x) x.^3 - (2+exp(1))*x.^2 + (2*exp(1)+1)*x + (1-exp(1)) - cosh(x-1);
x = linspace(0.5, 6.5, 100);
y = f(x);
figure(1);
plot(x,y)
title('f(x)=x.^3 - (2+exp(1))*x.^2 + (2*exp(1)+1)*x + (1-exp(1)) - cosh(x-1)');
xlabel('x');
ylabel('y');
grid on
y0 = zeros(100,1);
hold on
plot(x,y0)

%% Punto 2

disp('la prima radice appartiene all''intervallo [0.5, 1.5]');
if (f(0.5)*f(1.5)) < 0
    disp('si puo'' applicare il metodo di bisezione alla prima radice');
else
    disp('non si puo'' applicare il metodo di bisezione alla prima radice');
end

disp('la seconda radice appartiene all''intervallo [3, 4]');
if (f(3)*f(4)) < 0
    disp('si puo'' applicare il metodo di bisezione alla seconda radice');
else
    disp('non si puo'' applicare il metodo di bisezione alla seconda radice');
end

disp('la terza radice appartiene all''intervallo [6, 6.5]');
if (f(6)*f(6.5)) < 0
    disp('si puo'' applicare il metodo di bisezione alla terza radice');
else
    disp('non si puo'' applicare il metodo di bisezione alla terza radice');
end

%% Punto 3
% vedi funzione bisez.m

%% Punto 4

toll = 1e-12;

disp('calcolo seconda radice')
a1 = 3; b1 = 4;
[x1,it1]=bisez(a1,b1,toll,f);

disp('calcolo terza radice')
a2 = 6; b2 = 6.5;
[x2,it2]=bisez(a2,b2,toll,f);
