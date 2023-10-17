function dfun = derivata(fun,n)

% Calcola la derivata n-esima della funzione anonima scalare fun

if nargin == 1
	n = 1;
end

syms f(x)
f(x) = fun(x);
dfun = diff(f,n);
dfun = matlabFunction(dfun);
end