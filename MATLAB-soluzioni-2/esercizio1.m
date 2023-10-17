% Esercizio sui metodi di quadratura

% Punto 1

clear all
close all

figure(1);
a = 0; b = 0.8;
x_dis = linspace(a,b,1000);
f = @(x) cosh(x-0.5);
f_int = @(x) pi * (cosh(x-0.5)).^2;
y_dis=f(x_dis);
plot(x_dis,y_dis, 'r', 'linewidth',3 )
% opzioni grafico
grid on
axis([a b, 0 1.5])
set(gca,'FontSize',20)
title('f(x)=cosh(x-0.5)')
xlabel('x')
ylabel('f(x)')


%% Punto 2,3,4 -> vedi funzioni allegate

%% Punto 5

N = [1:20];
V_ex = pi*((sinh(1)+sinh(3/5))/4+2/5);
V_PMC = zeros(1,N(end));
V_TC = zeros(1,N(end));
V_SC = zeros(1,N(end));
for i = N
    V_PMC(i) = pmedcomp(a,b,i,f_int);
    V_TC(i) = trapcomp(a,b,i,f_int);
    V_SC(i) = simpcomp(a,b,i,f_int);
end
figure(2);
plot(N,V_PMC,'*',N,V_TC,'o',N,V_SC,'d',N, V_ex * ones(1,20), 'linewidth', 3)
grid on
legend('Pto medio composito','Trapezio composito','Simpson composito','V', 1) 
% opzioni grafico
set(gca,'FontSize',20)
xlabel('Numero sottointervalli')
ylabel('Integrale approssimato')
title('Andamento dell` integrale approssimato')

%% Punto 6

err_PMC= abs(V_ex - V_PMC);
err_TC= abs(V_ex - V_TC);
err_SC=abs(V_ex - V_SC);
H = (b-a)./N;
figure(3);
loglog(H,err_PMC,'*',H,err_TC,'o',H,err_SC,'d',H,H.^2,H,H.^4, 'linewidth', 3)
legend('err-PMC','err-TC','err-SC','H^2','H^4', 2)
% opzioni grafico
set(gca,'FontSize',20)
xlabel('Ampiezza dei sottointervalli')
ylabel('Errore')
title('Errore di quadratura in scala logaritmica')

rapp = err_PMC ./err_TC;
disp(rapp(end-4:end))

%% Punto 7

toll=1e-5;
d2f_int = @(x) 2*pi*(2*(cosh(x - 0.5)).^2 -1);
d4f_int = @(x) 8*pi*(2*(cosh(x - 0.5)).^2 -1);
K2 = max(abs(d2f_int(x_dis)));
K4 = max(abs(d4f_int(x_dis)));
Npmc = ceil(sqrt( ((b-a)^3 *K2)/(24*toll) ) )
Ntc = ceil(sqrt( ((b-a)^3*K2)/(12*toll) ) )
Nsc = ceil( ( ((b-a)^5*K4)/(16*180*toll) )^(1/4) )

VPM = pmedcomp(a,b,Npmc,f_int);
errPMCvero = abs(V_ex - VPM)
VTC = trapcomp(a,b,Ntc,f_int);
errTCvero = abs(V_ex - VTC)
VSC = simpcomp(a,b,Nsc,f_int);
errSCvero = abs(V_ex -VSC)

%% Punto 8 -> guardare la funzione gausscomp.m

%% Punto 9

V_GC = zeros(1,N(end));
for i = N
    V_GC(i) = gausscomp(a,b,i,f_int);
end
figure(4);
plot(N,V_GC,'o',N, V_ex * ones(1,20), 'linewidth', 3)
legend('Gauss composito','V', 3)
% opzioni grafico
set(gca,'FontSize',20)
xlabel('Numero sottointervalli')
ylabel('Integrale appr.')
title('Andamento dell` integrale approssimato')

err_GC= abs(V_ex - V_GC);
figure(5);
loglog(H,err_TC,'*',H,err_GC,'o',H,err_SC,'d',H,H.^2,H,H.^4, 'linewidth', 3)
legend('err-TC','err-GC','err-SC','H^2','H^4','Location','NorthWest')
% opzioni grafico
set(gca,'FontSize',20)
xlabel('Ampiezza dei sottointervalli')
ylabel('Errore')
title('Errore di quadratura in scala logaritmica')
