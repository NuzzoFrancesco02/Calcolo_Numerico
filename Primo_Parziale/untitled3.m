%{
x = linspace ( -3 ,3 ,1000); 
y = exp(-x.^2);
a = input('Digitare il valore di a:')




z = exp(-(x-a).^2);
if a>0
    s='Traslazione a destra';
elseif a<0
    s='Traslazione a sinistra';
else
    s='Nessuna traslazione';
end
k = input('Digitare il valore di k:')
w = exp(-(x.*k).^2);


if k>1
    t='Compressione orizzontale'; 
elseif k==1
    t='Nessun riscalamento'; 
elseif k>0 & k<1
    t='Dilatazione orizzontale'; 
else
    error('k deve essere maggiore di 0!'); 
end
%arrivano i grafici 
figure('name','Operazioni con le funzioni') 
subplot(3,1,1)
plot(x,y,'k') 
title('e^{-x^2}') 
grid on

subplot(3,1,2)
plot(x,z,'r')
grid on
title([s,' a= ',num2str(a)]) 
subplot(3,1,3)
plot(x,w,'b')
grid on
title([t,' k= ',num2str(k)])
%}
