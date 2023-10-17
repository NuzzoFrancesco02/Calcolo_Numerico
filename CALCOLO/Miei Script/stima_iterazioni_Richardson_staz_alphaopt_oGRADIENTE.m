
% STIMA VALIDA SOLO PER A SDP!!!!!!!!!!!

A = [];
P = eye(size(A));
FATTORE = 1000;		% FATTORE = 1/tol

if isequal(A,A')
	if eig(A) > 0
		fprintf('\nA è simmetrica definita positiva\n');
	else
		error('A è simmetrica ma non definita positiva');
	end
else
	error('A non è simmetrica');
end

K = spectralcond(P^-1*A);
d = (K - 1)/(K + 1);

%% Stima numero di iterazioni minimo per arrivare a tol
syms k;
Kmin = solve(d^k == 1/FATTORE,k)
Kmin = ceil(eval(Kmin))


%% Stima fattore abbattimento dopo k iteraazioni
k = 1;
tol = d^k