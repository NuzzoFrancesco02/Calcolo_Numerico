function[xk, k] = newtonvett(x0,fun,J,nmax,tol)

%Implementazione del metodo di newton per la ricerca degli zeri di una
%funzione vettoriale vicino ad x0

fprintf('\nMetodo di Newton per funzioni vettoriali:');

x = x0;
xk = [x];
err = tol + 1;
k = 0;

while(k<nmax && err>tol)
	Jx = J(x);
	if(round(det(Jx),12) == 0)
		error("Il metodo di Newton non è applicabile: Jacobiana con determinante == 0 all'iterazione %d", k);
	end
	k = k + 1;
	delta = - Jx \ fun(x);
	xn = x + delta;
	err = norm(xn - x);
	x = xn;
	xk = [xk, x];
end

if(k < nmax)
	fprintf("\nx_%d soddisfa la tolleranza sul residuo", k);
else
	fprintf("\nMax numero di iterazioni raggiunto! Errore %-6.4e", err);
end

fprintf("\nRadice calcolata :\n");
x
end