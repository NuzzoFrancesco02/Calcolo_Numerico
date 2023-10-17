function [t_h,u_h]=eulero_avanti(f,t_max,y_0,h)

% [t_h,u_h]=eulero_avanti(f,t_max,y_0,h)
%
% Risolve il problema di Cauchy 
%
% y'=f(y,t)
% y(0)=y_0
%
% utilizzando il metodo di Eulero in Avanti 
% (Eulero Esplicito) : u^(n+1)=u^n+h*f^n
%
% Input:
% -> f: function che descrive il problema di Cauchy 
%       (dichiarata ad esempio tramite inline o @). 
%       Deve ricevere in ingresso due argomenti: f=f(t,y)
% -> t_max: l'istante finale dell' intervallo temporale di soluzione 
%                 (l'istante iniziale e' t_0=0)
% -> y_0: il dato iniziale del problema di Cauchy
% -> delta_t: l'ampiezza del passo di discretizzazione temporale.
%
% Output:
% -> t_h: vettore degli istanti in cui si calcola la soluzione discreta
% -> u_h: la soluzione discreta calcolata nei nodi temporali t_h

% vettore degli istanti in cui risolvo la edo
t0=0;
t_h=t0:h:t_max;

% inizializzo il vettore che conterra' la soluzione discreta
N_istanti=length(t_h);
u_h=zeros(1,N_istanti);

% ciclo iterativo che calcola u_(n+1)=u_n+h*f_n
u_h(1)=y_0;

for it=2:N_istanti
    u_old=u_h(it-1);
    u_h(it)=u_old+h*f(t_h(it-1),u_old);
end

