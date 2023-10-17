%	Run dopo D = qrbasic(matrice)

for i = 1:length(D)-1
	v1(i) = abs(D1(i+1)/D1(i));
end

% "fattore di covergenza", più è basso meglio è
V1 = max(v1)