function[xk, k] = corde(x0,fun,a,b,nmax,tol)

%Implementazione del metodo delle corde per la ricerca degli zeri di una
%funzione vicino ad x0

fprintf('\nMetodo delle corde:');

x = x0;
xk = [x];
err = tol + 1;
k = 0;
q = (fun(b)-fun(a)) / (b-a);
while(k<nmax && err>tol)
	k = k + 1;
	x = x - fun(x)/q;
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