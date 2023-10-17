%% TEORIA
% e1*M1^p = e2*M2^p
% e1/h1^p = e2/h2^p
% K1*h1^p = K2*h2^p

%% ZERO-STABILITA'
% La zero-stabilità è una proprieta di un metodo numerico e permette di
% osservare il comportamento dello stesso in caso di perturbazioni numeriche 
% sulla soluzione in intervalli limitati. 
%
% Indichiamo con Nh il numero di intervalli in cui viene suddiviso
% l'intervallo [0,tf] con passo temporale h = tf/Nh, con g le dimensioni 
% del sistema, con z ∈ R_g la soluzione perturbata, con u ∈ R_g la soluzione 
% approssimata, con rho ∈ R_g un vettore che perturba tutte le soluzioni 
% all'istante di tempo n, per Eulero in avanti la soluzione approssimata è: 
% 
%   u(:,0) = y(:,0) + rho(:,0);
%   u(:,n+1) = u(:,n)+ h(f(t(n),u(:,n));
%
% la soluzione perturbata è: 
%
%   z(:,0) = y(:,0) + rho(:,0);
%   z(:,n+1) = z(:,n)+ h(f(t(n),z(:,n)+rho(:,n+1));
%
% Dove con la notazione (:,n) si vuole sottolineare la presenza di
% elementi vettoriali e non scalari.
% il metodo è zero-stabile se la differenza tra z(:,n) ed u(:,n) è controllata da una
% funzione della perturbazione vettoriale rho.
%
% In generale un metodo è zero-stabile se esistono h_0>0, C>0, eps_0>0, tali
% che, per ogni h ∈ (0,h_0] e per ogni eps ∈ (0,eps_0] 
%           ||rho(:,n)||<eps       per ogni n = 0,...,Nh
% implica:  max(||z(:,n)-u(:,n)||)<=C*eps      per ogni n = 0,...,Nh
% Dove C è una costante indipendente da h ma dipendente dalla lunghezza
% dell'intervallo |I| = tf - t0, e ||.|| è un'opportuna norma vettoriale.
%
% Nella forma A*y + g(t):
%       z(:,0) = y(:,0) + rho(:,0)
%       z(:,n+1) = x(:,n) + h*( A*z(:,n) + g(t(n)) + rho(:,n) )
%
%
% Nella forma A(y)*y + g(t):
%       z(:,0) = y(:,0) + rho(:,0)
%       z(:,n+1) = x(:,n) + h*( A(z(:,n))*z(:,n) + g(t(n)) + rho(:,n) )
%
%
% Nella forma A(t,y)*y + g(t):
%       z(:,0) = y(:,0) + rho(:,0)
%       z(:,n+1) = x(:,n) + h*( A(t(n),z(:,n))*z(:,n) + g(t(n)) + rho(:,n) )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORDINE DI CONVERGENZA
% Se l'errore associato ad un metodo numerico per EDO è tale per cui:
%               e(n)<=C*h^p, per n = 0,...,Nh
% allora:       e(n) ~ C*h^p --> log(e(n1)/e(n2)) ~ p * log_{h(1)/h(2)}(1) 
%               --> p ~= log(e(n1)/e(n2))/log(h(1)/h(2))
%
% Dalla teoria, poiché f(t,y) ∈ C2([t0,tf]) in y(i) con i = 1,...,ordine del sys
% Eulero in avanti converge con ordine p = 1 in h.

% Dalla teoria, poiché f(t,y) ∈ C2([t0,tf]) in y(i) con i = 1,...,ordine del sys
% Eulero indietro converge con ordine p = 1 in h.

% Dalla teoria, poiché f(t,y) ∈ C2([t0,tf]) in y(i) con i = 1,...,ordine del sys
% Crank-Nicolson converge con ordine p = 2 in h.

% Dalla teoria, poiché f(t,y) ∈ C2([t0,tf]) in y(i) con i = 1,...,ordine del sys
% Heun converge con ordine p = 2 in h.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASSOLUTA STABILITA'
% CASO COEFFICIENTI COSTANTI EULERO IN AVANTI
% In questo caso f(t,y) è differenziabile con continuità rispetto a 
% yi con i = 1,...,ord del sys e la jacobiana di f(t,y) coincide con la 
% matrice A, cerco l'autovalore massimo e il minimo di A: poiché sono entrambi
% finiti e a parte reale negativa, dalla teoria h_max è limitato 
% dall'autovalore massimo l_max. 
%
% Utilizzando la funzione di stabilità |R(z)|<1 con:
%                       z = 1 + h*l
% posso trovare h_max. In particolare per Eulero in avanti se l_max ∈ R 
% h_max = 2/abs(l_max), se l_max ∈ C -2*real(l_max)/abs(l_max)^2
%
%
% CASO COEFFICIENTI COSTANTI HEUN
% In questo caso f(t,y) è differenziabile con continuità rispetto a 
% yi con i = 1,...,ord del sys e la jacobiana di f(t,y) coincide con la 
% matrice A, cerco l'autovalore massimo e il minimo di A: poiché sono entrambi
% finiti e a parte reale negativa, dalla teoria h_max è limitato 
% dall'autovalore massimo l_max. 
%
% Utilizzando la funzione di stabilità |R(z)|<1 con:
%                       z = 1 + h*l + (h*l)^2/2
% posso trovare h_max
% CASO NON COSTANTI
% In questo caso f(t,y) è differenziabile con continuità rispetto a 
% yi con i = 1,...,ord del sys e la jacobiana di f(t,y) non coincide con la 
% matrice A. Pertanto, conoscendo dal punto precedente la soluzione esatta,
% costruisco la Jacobiana J(t) al variare di t da t0 a tf riempio due vettori: 
% v_max è il vettore che contiene gli autovalori massimi in ogni istante di tempo, 
% v_min gli autovalori minimi. Con l_max indico il massimo fra tutti gli istanti di
% tempo e l_min il minimo fra tutti gli istanti di tempo. Poiché sono entrambi
% finiti e a parte reale negativa, dalla teoria h_max è limitato dall'autovalore 
% massimo l_max. 
%
% Utilizzando la funzione di stabilità |R(z)|<1 con:
%                       EA: z = 1 + h*l
%                       HN: z = 1 + h*l + (h*l)^2/2
%% Numero condizionamento matrice di rigidezza

%Conoscendo gli autovalori della matrice è possibile calcolarne il numero
%di condizionamento:

%K(A) = C/h^2

%Notiamo quindi che diminuendo h, ovvero il passo aumenta la matrice di
%rigidezza e quindi la grandezza del sistema lineare, ma non solo il numero
%di condizionamento sarà sempre più elevato e dunque peggiore la soluzione
%del problema.

%Allo stesso tempo però l'errore per un problema ai limiti di
%Poisson-Dirichlet converge alla soluzione come una C*h^2. Questo risultato
% è dovuto all'approssimazione delle derivate tramite differenze finite di
% ordine 2.


%% Problema diffusione-trasporto-reazione

% Nel caso in cui eta>>mu il termine di trasporto è dominante rispetto alla
% diffusione (sto dando lo stesso peso all'informazione a monte e quella
% a valle, ma questo non è del tutto vero quando eta>>mu) 
% e questo comporta un problema, infatti per studiare correttamente tale
% problema con ordine 2 vi è imposto un bound sul passo h.

%  h < 2*mu/abs(eta)

% In caso non venga rispettato tale bound si può incorrere in oscillazioni
% numeriche sulla soluzione.


%% Upwind e Numero di Peclet

% Introducendo il numero adimensionale Peclet_locale = abs(eta)*h/(2*mu) e
% calcolando la soluzione approssimata ai nodi di u con ordine 2 notiamo
% che se Peh<1 si converge alla soluzione, mentre se Peh<1 numeratore e
% denominatore hanno segno discorde in modo alternato e questo causa le
% oscillazioni numeriche.

% Il problema è risolvibile rinunciando alla convergenza di ordine 2,
% approssimando il termine di trasporto u' con le differenze finite di
% primo ordine, il che porta alla formazione della viscosità artificiale
% ovvero un termine che è più grande della diffusione.