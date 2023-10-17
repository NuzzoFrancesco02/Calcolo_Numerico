function [t_h,u_h,iter_pf]=Crank_Nicolson(f,t_max,y_0,h)

%[t_h,u_h,iter_pf]=Crank_Nicolson(f,t_max,y_0,h)
%
% Risolve il problema di Cauchy 
%
% y'=f(y,t)
% y(0)=y_0
%
% utilizzando il metodo di Crank-Nicolson : u^{n+1}=u^n+1/2*h*( f^n +f^{n+1}). 
% Per ciascun istante temporale l'equazione per u^(n+1)(a priori non lineare) viene
% risolta con un metodo di punto fisso: u=u^n+1/2*h*( f^n + f(t^{n+1},u) ), 
% cioe' tramite la funzione di iterazione 
% 
% phi(x)=u^n+1/2*h*( f^n + f(t^{n+1},x) )
%
% la condizione per la convergenza del metodo di punto fisso ( |phi(x)'|<1) diventa:
%
% h*|\partial_x f(t^{n+1},x)| < 1
%
% ed e' garantita se h e' sufficientemente piccolo. 
%
% Input:
% -> f: function che descrive il problema di Cauchy (dichiarata ad esempio tramite 
%           inline o @) deve ricevere in ingresso due argomenti: f=f(t,y)
% -> t_max: l'istante finale dell' intervallo temporale di soluzione 
%                 (l'istante iniziale e' t_0=0)
% -> y_0: il dato iniziale del problema di cauchy: y(t=0)=dato_iniziale
% -> h: l'ampiezza del passo di discretizzazione temporale.
%
% Output:
% -> t_h: vettore contenente gli istanti in cui si calcola la soluzione discreta
% -> u_h: la soluzione discreta calcolata nei nodi temporali t_h
% -> iter_pf: vettore che contiene il numero di iterazioni di punto fisso utilizzate 
%                       per risolvere l'equazione non lineare ad ogni istante temporale.



% vettore degli istanti in cui risolvo la edo

t0=0;

t_h=t0:h:t_max;

% inizializzo il vettore che conterra' la soluzione discreta

N_istanti=length(t_h);

u_h=zeros(1,N_istanti);

% ciclo iterativo che calcola u^{n+1}=u^n+1/2*h*( f^n +f^{n+1}) . 
% Ad ogni iterazione temporale devo eseguire delle sottoiterazioni di punto 
% fisso per il calcolo di u^{n+1}: 
%
% u_(n+1)^(k+1) = u_n + h/2*( f^n + f_( t_(n+1) , u_(n+1)^k ) ). 

u_h(1)=y_0;

% parametri per le iterazioni di punto fisso
N_max=100;
toll=1e-5;

iter_pf=zeros(1,N_istanti);



for it=2:N_istanti
    
    % preparo le variabili per le sottoiterazioni
    
    u_old=u_h(it-1);
    f_old=f(t_h(it-1),u_old);

    t_pf=t_h(it);
        
    phi=@(u) u_old + 0.5*h * ( f_old + f( t_pf, u ) ) ;
    
    % sottoiterazioni
    
    [u_pf, it_pf] = ptofis(u_old, phi, N_max, toll);
    
    u_h(it)=u_pf(end);
    
    % tengo traccia dei valori di it_pf per valutare la convergenza delle 
    % iterazioni di punto fisso
 
    iter_pf(it)=it_pf;
    
end

