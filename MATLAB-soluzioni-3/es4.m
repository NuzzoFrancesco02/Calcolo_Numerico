%%%%%%%
% es4: PUNTO FISSO
%%%%%%%

clc; clear; close all

%% punto 1
f = @(x) cos(2*x).^2 - x.^2;
x = linspace(-pi/2,pi/2,1000);
plot(x,f(x),x,zeros(size(x)),'k');
grid on

%% punto 2 teorico. 
% phi' = 1 + Af'
% -1 < phi'(a) < 1  --> -2 < Af' < 0
% Notando che f'(a) < 0, il metodo di punto fisso 
% converge per 0 < A < -2/f'(a)

%% punto 3 - 1) determinazione dell'intervallo di ammissibilita' per A
figure
phi = @(x) x + 0.1*(cos(2*x).^2 - x.^2)
[succ1,it1] = ptofis(0.1,phi,1000,1e-10,-pi/2,pi/2);

df = @(x) -4*cos(2*x).*sin(2*x)-2*x
df_a = df(succ1(end)) % derivata di f nella soluzione = -2.7955
Asup = -2/(df_a)

%% punto 3 - 2) studio della convergenza al variare di A
% prove "random"
phi6 = @(x) x + 0.6*(cos(2*x).^2 - x.^2)
figure
[succ2,it2] = ptofis(0.1,phi6,1000,1e-10,-pi/2, pi/2);
% x0 lontano dalla soluzione
figure
[succ3,it3] = ptofis(2,phi6,1000,1e-10,-pi/2, pi);

% A > Asup
figure
phi75 = @(x) x + 0.75*(cos(2*x).^2 - x.^2)
[succ4,it4] = ptofis(0.1,phi75,1000,1e-10,-pi/2, pi);

%% punto 4
% ordine 1
% Il fattore di riduzione teorico e' phi'(a), ovvero 1 + A*(df_a)
[p,c] = stimap(succ1);
fprintf(' -----Stima teorica: %f\n\n',1 + 0.1*df_a)
[p,c] = stimap(succ2);
fprintf(' -----Stima teorica: %f\n\n',1 + 0.6*df_a)
% Perfetta aderenza tra il fattore di riduzione stimato e teorico 

%% punto 5
% Ho secondo ordine se phi'(a) = 0, ovvero se A = -1/f'(a)
A_opt = -1/df_a
phi_opt = @(x) x + A_opt*(cos(2*x).^2 - x.^2)
[succ5,it5] = ptofis(0.1,phi_opt,1000,1e-10,-pi/2, pi/2);
[p,c] = stimap(succ5); % ordine 2

%% punto 6
phiN = @(x) x - (cos(2*x).^2 - x.^2)./(-4*cos(2*x).*sin(2*x)-2*x)
figure
plot (x, phiN (x), 'b', [-3:0.01:3], [-3:0.01:3], 'g');
grid on
axis ([-2 2 -2 2]);

vett_sol=[];
vett_x0 = [0.01:0.01:1];
figure
for x0 = vett_x0
    [succ,it] = ptofis(x0,phiN,1000,1e-10,-pi/2, pi/2);
    vett_sol = [vett_sol succ(end)];
end
figure
plot(vett_x0,vett_sol,'o--')
