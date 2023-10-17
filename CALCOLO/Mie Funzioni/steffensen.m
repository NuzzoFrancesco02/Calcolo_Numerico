function[xk, k] = steffensen(x0,fun,nmax,tol)

%Implementazione del metodo di steffensen per la ricerca degli zeri di una
%funzione vicino ad x0

fprintf('\nMetodo di Steffensen:');

x = x0;
xk = [x];
err = tol + 1;
k = 0;

while(k<nmax && err>tol)
	k = k + 1;
	fx = fun(x);
	x = x - fx.^2./(fun(x+fx)-fx);
	err = abs(fun(x));
	xk = [xk; x];
end

if(k < nmax)
	fprintf("\nx_%d soddisfa la tolleranza sul residuo", k);
else
	fprintf("\nMax numero di iterazioni raggiunto! Errore %-6.4e", err);
end

fprintf("\nRadice calcolata = %-12.8f\n", x);
end