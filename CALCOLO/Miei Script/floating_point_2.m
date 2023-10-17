clear
clc

% SOLO PER BASE 2
b = 2;
m = [1 0 1 1];
e = [2];
s = [0];
t = length(m);

for i = 1:length(m)
	x(i) = ((-1)^s)*(b^e)*m(i)*b^(-i);
end
x = sum(x) + ((-1)^s)*(b^e)

%% eps da lunghezza mantissa
epsm = b^(-t)

%% t da eps
epsm = [];
t = -log2(epsm)

%% distanza tra un numero e il successivo
d = b^(1+e-t-1)