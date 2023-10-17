function [t_h,u_h,vett_it_pf]=eulero_indietro_pto_fisso(f,t_max,valore_iniziale,delta_t)
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
% delta_t           passo di integrazione
% 
%
% Parametri di uscita:
%
% t_h               vettore dei tempi in cui la soluzione viene calcolata
% u_h               vettore delle soluzioni calcolate in t_h
% vett_it_pf        vettore delle iterazioni del metodo di punto fisso 
%                   ad ogni passo

t0=0;

t_h=t0:delta_t:t_max;

% inizializzo il vettore che conterra' la soluzione discreta

N_istanti=length(t_h);

u_h=zeros(1,N_istanti);

% ciclo iterativo che calcola u_(n+1)=u_n+h*f_(n+1) . Ad ogni iterazione temporale devo eseguire delle sottoiterazioni
% di punto fisso per il calcolo di u_(n+1): u_(n+1)^(k+1) = u_n + h * f_( t_(n+1) , u_(n+1)^k ). Per garantire la
% convergenza di tale metodo la derivata della funzione di iterazione phi(x) u_n + h * f_( t_(n+1) , x ) deve avere
% derivata minore di 1 in modulo. Questo introdurra' delle condizioni su h.

u_h(1)=valore_iniziale;

% parametri per le iterazioni di punto fisso
N_max=100;
toll=1e-5;

vett_it_pf=zeros(1,N_istanti);



for it=2:N_istanti
    
    % preparo le variabili per le sottoiterazioni
    
    u_old=u_h(it-1);
    t_pf=t_h(it);
        
    phi=@(u) u_old + delta_t * f( t_pf, u );
    
    % sottoiterazioni
    
    [u_pf, it_pf] = ptofis(u_old, phi, N_max, toll);
    
    u_h(it)=u_pf(end);
    
    % tengo traccia dei valori di it_pf per valutare la convergenza delle iterazioni di punto fisso
    vett_it_pf(it)=it_pf;
    
end

