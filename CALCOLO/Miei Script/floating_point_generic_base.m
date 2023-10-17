clear
clc

b = 2;
m = [1 0 1 1];
e = [2];
s = [0];
t = length(m);

for i = 1:length(m)
	x(i) = ((-1)^s)*(b^e)*m(i)*b^(-i);
end
x = sum(x)

%% eps da lunghezza mantissa
epsm = b^(1-t)

%% t da eps
epsm = [];
t = - log(epsm)/log(b) + 1

%% distanza tra un numero e il successivo
d = b^(1+e-t)