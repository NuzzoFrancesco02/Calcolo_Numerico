clear; close all; clc;

%creo due autovettori ortonormali, che uso per creare la A sdp
v1=([4 3]/5)';
v2=([-3 4]/5)';

%questi sono i due autovalori associati
l1=5;
l2=10;

V  = [v1,v2]; 
D1 = diag([l1,l2]);
A1 = V*D1*V';

%se cambio gli autovalori la forma quadratica associata e' piu' schiacciata
%(NB autoval ~ 1/semiasse percio' autoval piccolo =>  semiasse lungo ), come
%si vede qualche riga sotto, quando disegniamo le forme quadratiche Phi1 e
%Phi2
l1=1;
l2=10;

D2 = diag([l1,l2]);
A2 = V*D2*V';

%questo e' il termine noto della forma quadratica che vogliamo minimizzare
b=[4 8]';

%definisco la griglia 
x=linspace(-10,10,80);
y=linspace(-10,10,80);

[X,Y]=meshgrid(x,y);

%definisco le forme quadratiche Phi1 e Phi2
Phi1 =  0.5* ( A1(1,1)*X.^2 + A1(2,2)*Y.^2 + 2*A1(1,2)*X.*Y ) - b(1)*X - b(2)*Y;
Phi2 =  0.5* ( A2(1,1)*X.^2 + A2(2,2)*Y.^2 + 2*A2(1,2)*X.*Y ) - b(1)*X - b(2)*Y;

% Avrei potuto definire le forme quadratiche anche tramite anonymous
% function:
% Phi1 =  @(X,Y) 0.5* ( A1(1,1)*X.^2 + A1(2,2)*Y.^2 + 2*A1(1,2)*X.*Y ) - b(1)*X - b(2)*Y;
% Phi2 =  @(X,Y) 0.5* ( A2(1,1)*X.^2 + A2(2,2)*Y.^2 + 2*A2(1,2)*X.*Y ) - b(1)*X - b(2)*Y;
% Per poi valutarle nelle matrici X ed Y:
% Phi1_val = Phi1(X, Y)
% Phi2_val = Phi2(X, Y)
% ed infine usare nei comando surf e contour Phi1_val e Phi2_val.


%confrontiamo graficamente le due forme quadratiche 
%1) questo crea il grafico 3D
surf(X,Y,Phi1,'Lines','no');
title('phi 1')

hold on
%2) aggiungo al grafico 3D le linee di livello; l'ultimo parametro da' il numero di linee di livello
%desiderate
contour(X,Y,Phi1,20)

%rifaccio lo stesso per la seconda forma quadratica
figure;
surf(X,Y,Phi2,'Lines','no');
title('phi 2')
hold on;
contour(X,Y,Phi2,20)

%per rendersi bene conto del fatto che la seconda forma quadratica e' piu'
%schiacciata la cosa migliore e' confrontare le linee di livello delle due
%forme quadratiche

figure
contour(X,Y,Phi1,20)
title('phi 1 - linee di livello')

figure
contour(X,Y,Phi2,20)
title('phi 2 - linee di livello')


%% risolvo il problema di minimizzazione della forma phi2 con il metodo di Richardson stazionario (alpha scelto in modo da convergere)
x0=[-9 -9]';
nmax=1000;
toll=1e-7;
[Sg_rich,n_it_rich] = richardson_it(A2,b,eye(2),x0,toll,nmax,0.05);


%% plotto la storia di convergenza
figure;
contour(X,Y,Phi2,80)
title('\alpha = 0.05')
%con questo comando il grafico avra' la stessa scala per l'asse x e per l'asse y
axis equal
hold on;
plot(Sg_rich(1,:),Sg_rich(2,:),'-or','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r')


%% risolvo il problema di minimizzazione della forma phi2 con il metodo di Richardson stazionario (alpha scelto in modo da divergere)
x0=[-9 -9]';
nmax=1000;
toll=1e-7;
[Sg_rich2,n_it_rich2] = richardson_it(A2,b,eye(2),x0,toll,nmax,0.24);


%% plotto la storia di convergenza
figure;
contour(X,Y,Phi2,80)
title('\alpha = 0.24')
axis equal
hold on;
plot(Sg_rich2(1,1:6),Sg_rich2(2,1:6),'-or','LineWidth',2)



%% risolvo il problema di minimizzazione della forma quadratica phi2 con il metodo del gradiente
x0=[-9 -9]';
nmax=100;
toll=1e-7;
[Sg2,n_it2] = richardson_it(A2,b,eye(2),x0,toll,nmax);


%% plotto la storia di convergenza
figure;
contour(X,Y,Phi2,80)
title('Gradiente')
axis equal
hold on;
plot(Sg2(1,:),Sg2(2,:),'-or','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r')


%% voglio risolvere il problema di minimizzazione della forma quadratica phi2 con il
%  metodo del gradiente precondizionato, tramite il precondizionatore P

x0=[-9 -9]';
nmax=100;
toll=1e-7;

%P e' una matrice che approssima A2: al posto di avere l2=10 usiamo l2=5
P=V *diag([l1,5])*V';

%quello che vogliamo ottenere quindi e' una nuova forma quadratica Phiprec meno
%schiacciata, ma che condivida il minimo con Phi2; 

%in effetti la forma quadratica precondizionata e' meno schiacciata, ma ha
%lo stesso minimo

Aprec=P\A2;
eig(Aprec)
bprec=P\b;

Phiprec = 0.5* ( Aprec(1,1)*X.^2 + Aprec(2,2)*Y.^2 + 2*Aprec(1,2)*X.*Y ) ...
    - bprec(1)*X - bprec(2)*Y;

figure
contour(X,Y,Phi2,80)
grid on
title('phi 2 - linee di livello')

figure
contour(X,Y,Phiprec,80)
grid on
title('phi prec - linee di livello')

%% risolvo con il gradiente precondizionato
[Sg_prec,n_it_prec] = richardson_it(A2,b,P,x0,toll,nmax);



%% plotto la storia di convergenza
figure;
contour(X,Y,Phi2,60)
axis equal
hold on;
title('Gradiente Precondiz. vs Gradiente')

%disegno in grigio le direzioni del gradiente normale
plot(Sg2(1,:),Sg2(2,:),'--','LineWidth',2,'Color',[0.48 0.48 0.48])
 
%poi sovrappongo in rosso quelle del gradiente precondizionato
plot(Sg_prec(1,:),Sg_prec(2,:),'-or','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r')



%% 
% applicare il gradiente precondizionato e' equivalente ad applicare il
% gradiente non precondizionate alla Phiprec? No! 
[Sg_phiprec,n_it_phiprec] = richardson_it(Aprec,bprec,eye(2),x0,toll,nmax);

figure;
contour(X,Y,Phiprec,60)
axis equal
title('Gradiente su P^{-1}A')
hold on;
plot(Sg_phiprec(1,:),Sg_phiprec(2,:),'-or','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r')


%% verifiche di ortogonalita'

%direzioni del gradiente precondizionato
dir_prec=Sg_prec(:,2:5)-Sg_prec(:,1:4)

%direzioni del gradiente applicato al sistema precondizionato
dir_phiprec=Sg_phiprec(:,2:5)-Sg_phiprec(:,1:4)

%le direzioni del gradiente applicato al sistema precondizionato devono essere ortogonali
disp('gradiente applicato al sistema precondizionato: test ortogonalita'' incrementi consecutivi (se 0 allora ortogonali, altrimenti non ortogonali)')
dir_phiprec(:,1)'*dir_phiprec(:,2)
dir_phiprec(:,2)'*dir_phiprec(:,3)
dir_phiprec(:,3)'*dir_phiprec(:,4)

%mentre quelle del gradiente precondizionato in generale non lo saranno ...
disp('gradiente precondizionato: test ortogonalita'' incrementi consecutivi (se 0 allora ortogonali, altrimenti non ortogonali)')
dir_prec(:,1)'*dir_prec(:,2)
dir_prec(:,2)'*dir_prec(:,3)
dir_prec(:,3)'*dir_prec(:,4)

%... ma sono a due a due P ortogonali:
disp('gradiente precondizionato: test P-ortogonalita'' incrementi consecutivi (se 0 allora ortogonali, altrimenti non ortogonali)')
dir_prec(:,1)'*P*dir_prec(:,2)
dir_prec(:,2)'*P*dir_prec(:,3)
dir_prec(:,3)'*P*dir_prec(:,4)

%direzioni non consecutive invece in generale non saranno P ortogonali  
disp('gradiente precondizionato: test P-ortogonalita'' incrementi non consecutivi (se 0 allora ortogonali, altrimenti non ortogonali)')
dir_prec(:,1)'*P*dir_prec(:,3)
dir_prec(:,1)'*P*dir_prec(:,4)
dir_prec(:,2)'*P*dir_prec(:,4)




%% possiamo fare di meglio, ed avere tutte direzioni A2-ortogonali?? usiamo il gradiente coniugato
[Sg_conj,n_it_conj] = conjgrad_it(A2,b,x0,nmax,toll);

figure;
contour(X,Y,Phi2,60)
axis equal
hold on;
title('Gradiente Coniugato vs Precondizionato')

%plotto in grigio le iter di gradiente precondizionato
plot(Sg_prec(1,:),Sg_prec(2,:),'--','LineWidth',2,'Color',[0.48 0.48 0.48])
plot(Sg_conj(1,:),Sg_conj(2,:),'-or','LineWidth',2)
