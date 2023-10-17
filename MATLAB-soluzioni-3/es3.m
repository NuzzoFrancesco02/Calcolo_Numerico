clear
clc
close all

%% Punto 1

a = -1;  b = 6;
x = linspace( a, b, 1000 );
f = @(x) atan(7*( x - pi/2)) + sin((x-pi/2).^3);
y = f(x);
y0 = zeros(length(x),1);
figure(1);
plot(x,y,'b',x,y0,'r')
title('f(x) = atan(7*(x - pi/2)) + sin((x-pi/2)^3)')
xlabel('x')
ylabel('y')
legend('y = f(x)', 'y = 0', 'Location', 'SouthEast')
grid on;

%% Punto 2

alpha = pi/2; % valore esatto dello zero di f(x)

df = @(x) 7 ./ ( 1 + 49 * ( x-pi/2 ).^2 ) + 3 * (x-pi/2).^2 .* cos( (x-pi/2).^3 );

% controllo che alpha = pi/2 sia uno zero semplice
if (df(alpha) ~= 0)
    disp('lo zero alpha = pi/2  e'' semplice')
end

nmax = 1000;
toll = 1e-10;

disp('Metodo di Newton con X0 = 1.5')
x0 = 1.5;
[xvect,it] = newton(x0,nmax,toll,f,df);  
err = abs(xvect(end)-alpha)

disp('Metodo di Newton con X0 = 4')
x0 = 4;
[xvect,it] = newton(x0,nmax,toll,f,df);  
err = abs(xvect(end)-alpha)

%% Punto 3
if (f(a)*f(b)) < 0
    disp('si puo'' applicare il metodo di bisezione');
    nmax_b = 1000;
    toll_b = ( b - a ) / (2^31);
    [xvect,it] = bisez(a,b,toll_b,f);
    err = abs(xvect(end)-alpha)
else
    disp('non si puo'' applicare il metodo di bisezione');
end


%% Punto 4
% vedi biseznewton.m

%% Punto 5

nmax_b = 5;
nmax_n = 1000;
[xvect,it] = biseznewton(a,b,nmax_b,nmax_n,toll,f,df);

disp('l''errore finale e''')
disp(abs(xvect(end)-alpha))