function Z = allzeros(f,a,b,N)

% Trova tutti gli zeri di f (scalare) nell'intervallo [a,b]

if nargin < 4
	N = 100;
end

if length(a) == length(b)
	n = length(a);
else
	error('a e b hanno dimensioni discordi');
end

for i = 1:n
	temp = a(i);
	a(i) = min(a(i),b(i));
	b(i) = max(temp,b(i));
end

dx = (b-a)./N;
x2 = a;
y2 = f(x2);
Z = [];

for i = 1:N
	x1 = x2;
	y1 = y2;
	x2 = a + i*dx;
	y2 = f(x2);
	if y1*y2 <= 0				% Teorema di Rolle
		Z = [Z; fsolve(f,(x2*y1-x1*y2)/(y1-y2))];	% Approssimazione lineare per essere il più vicino possibile
	end
end
end