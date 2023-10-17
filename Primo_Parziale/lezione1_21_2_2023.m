clear
clc
%{
%Esercizio 1.1
v1 = 2.^(0:10);


v2 = cos(pi./(1:10))';


v3 = 0.1./(2.^(0:5));


% v4=zeros(1,19);
% j=1; k=6;
% for i = 1 : 19
%     if(mod(i,2) == 0)
%         v4(i) = 0;
%         fprintf("0\n");
%     else
%         v4(i) = exp(j)+(-1).^(j+1).*k;
%         fprintf("v[%d]=e^%d +(%d)*%d \n",i , j,((-1)^(j+1)),k);
%         j=j+1; k=k+5;
%     end
% end
v4 = zeros(1,19);
v4_dispari = exp(1:10)-(-1).^(1:10)+(6:5:51);
v4(1:2:20) = v4_dispari
%}







%{
Esercizio 1.2
A1 = M(1:3,1:3);
A2 = M((1:2:5),[1,2,4]);
A3 = M((2:4),[1,3,4]);
M = 2*eye(5) + diag(5*ones(4,1),-1) + diag(10*ones(3,1),-2) + diag(10*ones(3,1), 2) + diag(40,4) + diag(40,-4);
sum(M,"all");
%}



%{
%Esercizio 1.3
B = diag(ones(1,10)) + [0, ones(1,8), 0]' * [1, zeros(1,8), 1] + [1, zeros(1,8), 1]' * [0, ones(1,8), 0];


C = diag([1:200]) + diag(ones(199,1),1) + diag(ones(199,1),-1) + diag(0.5*ones(198,1),2) + diag(0.5*ones(198,1),-2);


D = diag([20:-2:2]) + diag(0.5*ones(1,8), -2) + diag(3.^[0:8],1); format("default");
%}


%{
%Esercizio 1.4
x = [0:3];
function1_eval = 'x.*sin(x)+(1/2).^nthroot(x,2)';
function1_inline = inline('x.*sin(x)+(1/2).^nthroot(x,2)', 'x');
function1_anon = @(x) x.*sin(x)+(1/2).^nthroot(x,2);
eval(function1_eval);
function1_inline(x);
function1_anon(x);

function2_eval = 'x.^4 + log(x.^3+1)';
function2_inline = inline('x.^4 + log(x.^3+1)','x');
function2_anon = @(x) x.^4 + log(x.^3+1);
eval(function2_eval);
function2_inline(x);
function2_anon(x);
%}


%{
%Esercizio 1.5
clear x;
x=[0:0.01:6];
f1 = @(x) 2+(x-3).*sin(5.*(x-3));
plot(x,f1(x), 'k');
hold on;
plot(x,-x+5,'--k');
hold on;
plot(x,x-1,'--k');


clear x;
hold off;
x=[0.1:0.01:10];
f2 = @(x) log(x).^2;
plot(x,f2(x),'k');
%}


%{
%Esercizio 4.1
A = mat_hilbert(5);
M = hilb(5);
if(A == M)
    fprintf('Correct!\n');
else
    error('Sbagliato!');
end
%}


%Esercizio 4.2
%succ = treuno(7);


%{
%Esercizio 4.3
anni = 0;
conto = 10000;
while(conto<(10^6))
    conto = conto + conto*0.02;
    conto = conto + 10000;
    anni = anni +1;
end
fprintf('\nOccorrono %d anni\n', anni);
%}




%{

%Esercizio 4.4
tic
n=100;
coppie = rand(n,2);
m=0.00;
for i = 1:n
    if(norm(coppie(i,1:2))<=1)
        m = m+1;
    end
end
pigreco = 4*m/n
toc
%}



%{
%Esercizio 4.5

sum=0; i=1;
while(sum < 88)
    sum = sum + i;
    i = i+1;
end
i-1
%}



%{
%Esercizio 4.6
A = zeros(10,1);
for k = 1:10
    if (k == 2 || k == 6)
        A(k,1) = 1/k;
    else A(k,1) = 1/(k-2)*(k-6);
    end
end
%}




%{
%Esercizio 4.7
v = zeros(1,10);
for k= 0:8
    v(1,k+1)=(2*k+1)^2;
end
%}



%{
%Esercizio 4.8

function vk=v_fun(n)
        vk = zeros(1,10);
   for k= 0:n
        vk(1,k+1)=(2*k+1)^2;
    end
end 
%}


%{
%Esercizio 4.9
r = erone(100);
x = 1:length(r);

plot(x,r,'-*r',x,sqrt(100).*ones(1,length(r)),'b');
function [r] = erone(n)
tol = 1e-6;
r(1) = n;
err = 1;
k = 1;
while (err > tol)
    k = k+1;
    r(k) = 0.5*(r(k-1) + n/r(k-1));
    err = abs(r(k)-r(k-1));
end
end

%}


%{
%Esercizio 4.10
p = zeros(1,30); q = p; b = p; a = p;
x = (1:1:length(p));
a(1) = sqrt(2); b(1) = 2; p(1) = 2*sqrt(2); q(1) = 4;
for i = 2: length(p)
    a(i) = sqrt(2).*(sqrt(1-sqrt(1-((a(i-1).^2)./4))));
    b(i) = a(i) ./(sqrt(1 -a(i-1) .^2 ./4));
    p(i) = (2^(i)).*a(i);
    q(i) = (2^(i)).*b(i);
end
plot(x,p,'-*r',x,q,'-*b');
hold on;
%}



%{
%Esercizio 4.12
x=(-10:0.01:10);

z = (-sqrt(x.^2-x)).*(x<0)+(-x.^2+2.*x).*exp(-x).*(x>=0);
subplot(1,2,1);
plot(x,z);
subplot(1,2,2);
plot(x,f4_12(x));

function f = f4_12(x)
for i = 1:length(x)
    if(x(i)<0)
        f(i) = -sqrt(x(i).^2-x(i));
    else
        f(i) = (-x(i).^2+2.*x(i)).*exp(-x(i));
    end
end
end
%}




%{
%Esercizio 4.13
M = matrice;
function M = matrice()
n = input('\nInserisci dimensione matrice: ');
M = zeros(n);
k = input('\nInserisci penultima cifra del codocie matricola: ');
for i = 1 : n
    for j = 1 : n
        M(i,j) = 2*i*j+(k+1);
    end
end
end
%}

%Esercizio 4.14
%f4_14 = @(x,y) atan(y./x)-(sin(x.*sqrt(abs(y))).^2);



%{
%Esercizio 4.14
[x,y]= meshgrid(0.1:0.1:10);
z = atan(y./x)-(sin(x.*sqrt(abs(y))).^2);
surfl(x, y, z);
shading faceted;
%}


%{
%Esercizio 4.15
r=2.5+mod(1,2);
succ4_15=zeros(1,100);
succ4_15(1,1)=0.5;
band=0;
for i = 1:99
    succ4_15(1,i+1) = r.*succ4_15(1,i).*(1-succ4_15(1,i));
    if (abs(succ4_15(1,i+1)-succ4_15(1,i))<(1.e-3) && band == 0)
       band=1; 
       num=succ4_15(1,i+1);
       n=i+1;
    end
end
if(band==1)
    fprintf('\na_%d = %f è il primo numero che soddisfa la condizione\n',n,num);
else
    fprintf('\nNon esiste nessun numero che soddisfa la condizione');
end
%}



%{
%Esercizio 4.17
%Bubble sort
vett4_17=rand(1,20);
for i = 1:(length(vett4_17)+1)
    for j = 1:length(vett4_17)-i
        if(vett4_17(j)<vett4_17(j+1))
            temp=vett4_17(j);
            vett4_17(j)=vett4_17(j+1);
            vett4_17(j+1)=temp;
        end
    end
end
%}


%{
%Esercizio 4.19

x=(-pi:0.1:pi);
plot(x,f4_19(x));
function f = f4_19(x)
for i = 1: length(x)
    if ((x(i)>=-pi) && (x(i)<0))
        f(i)=sin(x(i));
    end
    if ((x(i)>=0) && (x(i)<=pi))
            f(i)=(x(i).^2+x(i))./6;
    end
end
end

%}



%Esercizio 4.20
%{
n = input('\nInserisci un numero: ');
if((floor(sqrt(n))-sqrt(n))==0)
    fprintf('\n %d è un quadrato perfetto',n);
else fprintf('\n %d non è un quadrato perfetto\n', n);
end
%}



%{
%Esercizio 4.23
Prova = matrix_diagonals()
function M=matrix_diagonals()
    n = input('\nInserisci lato matrice: ');
    M=ones(n);
    for i= 1:2:n-1
        M=M+diag(-1.*ones(1,n-i),i);
        M=M+diag(-1.*ones(1,n-i),-i);
    end
end 
%}



%{
%Esercizio 4.26
n = zeros(1,30);
n(1) = 1;
for i = 2:length(n)
    n(i) = ((n(i-1).^2)+2)./(2.*n(i-1));
end
x = [1:length(n)];
plot(x,n,'*-b');
hold on;
sqr=sqrt(2).*ones(length(x));
plot(x,sqr,'r');
%}



%{
%Esercizio 4.27

A = npern()
function M = npern()
    n = input('\nInserisci dimensione matrice: ');
    M = zeros(n);
    for i = 1:n
        M = M+diag(i.*ones(1,n-i+1),i-1);
    end
end

%}

%{
%Esercizio 4.27
fx=@(x) (x.^2+3)./(x-1);
x=[-100:0.01:100];
plot(x,fx(x),'r');
grid on;
axis([-100, 100, -100, 100]);
%}


%{
x = -5:0.01:5;
f = @(x) sqrt(9.*(((-x.^2)./4)+1));
plot (x,f(x),x,-f(x));


f1 = @(x) (5.*ones(1,length(x)));
plot(x,f1(x),x,-f1(x),f1(x),x,-f1(x),x);
%}

%{

k = input('\nInserisci ordine di approssimazione: ');
tic;
x = (-2*pi : 0.1 : 2*pi);
esp = 1;
f = zeros(1,length(x));
for i = 1 :2: k
    f = f + ((-1)^(esp-1).*((x.^i)/factorial(i)));
    esp = esp + 1;
    if(mod(i,2) ~= 0)
    %plot(x,f,'LineWidth',3,'r',x,sin(x),'LineWidth',0.5,'--k');
    plot(x,f,'Linewidth', 2, 'DisplayName', [num2str(i),'° ordine']);
    legend
    hold on;
    plot(x,sin(x), 'Linewidth', 0.5);
    grid on;
    axis([min(x), max(x), -5, 5]);
    end
end
tempo = toc;
fprintf('\nIl programma ha impiegato: %f secondi', tempo);
%plot(x,f,x,sin(x))
%}


%{

funzione = input('\nInserisci funzione: ', 's');
x0 = input('\nInserisci punto di passaggio: ');
k = input('\nInserisci ordine di approssimazione: ');
x = -10:0.01:10;
h = 0.01;
f = inline(funzione, 'x');
risultato = f(x0).*ones(1,length(x));
der = f(x);
%px_str = num2str(f(x0));
df = @(x) diff(x)./h;
j = find(x==x0);

for i = 1 : k
    der = df(der);
    px = @(t) der(j).*((t-x0).^i)./factorial(i);
    %px_str = strcat(px_str,[ '+' num2str(der(j)) '.*((t-' num2str(x0) ').^' num2str(i) ')./' num2str(factorial(i))]);
    risultato = risultato + px(x);
end

%px = inline(px_str, 't');
plot(x,f(x),x,risultato);
hold on;
axis([-5 5 -pi pi]);
%}

%

%}






 