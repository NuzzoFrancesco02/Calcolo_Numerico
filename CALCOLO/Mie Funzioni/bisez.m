function[xk,k] = bisez(a,b,fun,tol)

%Implementazione del metodo di bisezione per la ricerca degli zeri di una
%funzione in un intervallo [a,b]

%nmax = [log2((b-a)/tol)-1]

if(fun(a)*fun(b)>0)
	error("Il metodo di bisezione non è applicabile nell'intervallo [%d,%d] perchè assume lo stesso segno negli estremi", a, b);
end

xk = [];
err = tol + 1;
k = -1;
nmax = ceil(log2((b-a)/tol) - 1);

while(k<nmax && err>tol)
	k = k + 1;
	x = (a+b)/2;
	xk = [xk; x];
	fx = fun(x);
	if(fun(a)*fx < 0)
		b = x;
	else
		a = x;
	end
	
	err = abs(fx); 
end

if(k < nmax)
	fprintf("\nx_%d soddisfa la tolleranza sul residuo", k);
else
	fprintf("\nMax numero di iterazioni raggiunto! Errore %-6.4e", err);
end

fprintf("\nRadice calcolata = %-12.8f\n", x);
end
%		% ="variabile"		- ="allinea a sx"		 00 ="spazio da lasciare"
%		. = "non ricordo"	00 ="numero di cifre decimali"		x ="exp/float/etc..."