function intfun = integrale(fun)

% Calcola l'integrale indefinito della funzione anonima scalare fun

syms f(x)
f(x) = fun(x);
intfun = int(f);
intfun = matlabFunction(intfun);
end