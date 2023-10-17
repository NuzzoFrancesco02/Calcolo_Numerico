%{
DOMANDA 1
in ogni intervallo [2^e, 2^(e+1)] sono distanziati di 2^(e-t)
e = 1, valore successivo a 2;
t = 4;
2^(e-t) = 2^(1-4)

1 e poi 1+em, em = beta^(1-t) o beta^(-t) se normalizzato o no
distanza = successivo di 2 : 2*(1 + em)-2 = 2*(1 + beta^(1-t))-2

se normalizzato 

DOMANDA 2
s = 1
m = 10101 in base 2
e = 1

Convenzione A = x = (-1)^s * 2^e *somm(i=1 a s)(ai*2^-i) = -1.3125
Convenzione B = x = (-1)^s * 2^e * (1+ mtilde) = -2*(1+1/2+1/8*1/32) =
-3.7125
%}

%%
f1 = @(x) (x-1).^7;
f2 = @(x) x.^7 - 7.*x.^6 + 21.*x.^5 - 35.*x.^4 + 35.*x.^3 - 21.*x.^2 + 7.*x -1;
x = 1.01;
err = 100*(abs(f1(x)-f2(x))/abs(f1(x)));

%%
n = 100;
A = zeros(n);
for i = 1:n
    for j = 1:n
        k = i - j;
        A(i,j) = 100 - abs(k);
    end 
end
s1 = sum(A(100,:))
s2 = sum(diag(fliplr(A)))

%%
n = 10;
S = 1e5;
x = zeros(1,n+1);
x(1) = S;
for i = 2:n+1
    x(i) = (x(i-1) + S/x(i-1))/2;
end 
plot(1:n+1,x)
x(11)

