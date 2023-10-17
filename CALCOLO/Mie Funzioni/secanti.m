function[xk, k] = secanti(x0,x1,fun,nmax,tol)

%Implementazione del metodo delle secanti per la ricerca degli zeri di una
%funzione vicino ad x0

fprintf('\nMetodo delle secanti:');

xk = [x0; x1];
x = x1;
err = tol + 1;
k = 1;

while(k<nmax && err>tol)
	k = k + 1;
	fx = fun(x);
	qk = (fx - fun(xk(end-1)))/(xk(end) - xk(end-1));
	xn = x - fx/qk;
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