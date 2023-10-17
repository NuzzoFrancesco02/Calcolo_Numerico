function [t_h,u_h]=eulero_indietro_sistemi(f,t_max,valore_iniziale,h)
% [t_h,u_h,vett_it_newton]=eulero_indietro_pto_fisso (f, t_max, valore_iniziale, delta_t)
%
% Metodo di Eulero all'indietro per la soluzione di equazioni differenziali
% ordinarie con metodo di punto fisso per il calcolo della soluzione ad ogni
% step
%
% Parametri di ingresso:
%
% f                 funzione di t, u associata al problema di Cauchy
% t_max             estremo di integrazione
% valore_iniziale   valore della soluzione al tempo iniziale t0 = 0
% h         passo di integrazione
% 
%
% Parametri di uscita:
%
% t_h               vettore dei tempi in cui la soluzione viene calcolata
% u_h               vettore delle soluzioni calcolate in t_h
% vett_it_pf        vettore delle iterazioni del metodo di punto fisso 
%                   ad ogni passo

t0=0;

t_h=t0:h:t_max;

% inizializzo il vettore che conterra' la soluzione discreta

N_istanti=t_max/h;

u_h=zeros(length(valore_iniziale),N_istanti+1);

% ciclo iterativo che calcola u_(n+1)=u_n+h*f_(n+1) . Ad ogni iterazione temporale devo eseguire delle sottoiterazioni
% di punto fisso per il calcolo di u_(n+1): u_(n+1)^(k+1) = u_n + h * f_( t_(n+1) , u_(n+1)^k ). Per garantire la
% convergenza di tale metodo la derivata della funzione di iterazione phi(x) u_n + h * f_( t_(n+1) , x ) deve avere
% derivata minore di 1 in modulo. Questo introdurra' delle condizioni su h.

u_h(:,1)=valore_iniziale;

% parametri per le iterazioni di punto fisso
N_max=100;
toll=1e-5;


for it=2:N_istanti+1
    
    % preparo le variabili per le sottoiterazioni
    
    u_old=u_h(:,it-1);
    t_pf=t_h(it);
        
    phi=@(u) u_old + h* f( t_pf, u );
    
    % sottoiterazioni
    
    [u_pf] = ptofis_vettoriale(u_old, phi, N_max, toll);
    
    u_h(:,it)=u_pf(:,end);
    
    
 end
end

