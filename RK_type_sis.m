function [th,uh]=RK_type_sis(f,B,t0,tmax,y0,h) %B è l'array di Butcher
%Variante per sistemi

%Si assume che i Ki dipendano solo dai Ki già precedentemente calcolati
%perché non saprei come implementare un caso più generico (quindi K1
%dipende solo da t0 e da u0, K2 da t0,u0 e K1 e così via)

%IPOTESI 
%a_ij=0 per ogni j>=i (B è triangolare inferiore)

%IMPORTANTE 
%CONTROLLARE CHE L'ARRAY DI BUTCHER RISPETTI LE IPOTESI, PERCHÉ SONO SOLO
%TEORICHE AL FINE DI AVERE UN RISULTATO CORRETTO, MA LA FUNCTION FUNZIONA
%UGUALMENTE (restituendo un risultato sbagliato) 

th=t0:h:tmax;
N=length(th);

a=B(1:end-1,2:end);
b=B(end,2:end);
c=B(1:end-1,1);

K=zeros(size(y0,1),length(c));
prod1=zeros(size(y0,1),1);
prod2=zeros(size(y0,1),1);

uh=zeros(size(y0,1),N);
uh(:,1)=y0;

for k=2:N
    for i=1:length(c)
        for j=1:size(y0,1)
            prod1(j)=a(i,:)*K(j,:)';
        end
        K(:,i)=f(th(k-1)+h*c(i),uh(:,k-1)+h*prod1);
    end
    for n=1:size(y0,1)
        prod2(n)=b*K(n,:)';
    end
    uh(:,k)=uh(:,k-1)+h*prod2;
end
end