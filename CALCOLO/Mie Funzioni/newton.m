function[xk, k] = newton(x0,fun,dfun,nmax,tol,mol)

%Implementazione del metodo di newton per la ricerca degli zeri di una
%funzione vicino ad x0

if(nargin == 5)
	fprintf('\nMetodo di Newton:');
	mol = 1;
else
	fprintf('\nMetodo di Newton modificato:');
end

x = x0;
xk = [x];
err = tol + 1;
k = 0;

while(k<nmax && err>tol)
	dfx = dfun(x);
	if(dfx == 0)
		error("Il metodo di Newton non è applicabile: derivata nulla nel punto %d all'iterazione %d", x0, k);
	end
	k = k + 1;
	xn = x - mol*fun(x)/dfx;
	err = abs(xn - x);
	x = xn;
	xk = [xk; x];
end

if(k < nmax)
	fprintf("\nx_%d soddisfa la tolleranza sul residuo", k);
else
	fprintf("\nMax numero di iterazioni raggiunto! Errore %-6.4e", err);
end

fprintf("\nRadice calcolata = %-12.8f\n", x);
end