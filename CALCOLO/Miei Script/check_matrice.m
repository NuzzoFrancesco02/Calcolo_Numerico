A = [];

%% Singolare
if round(det(A),8)
    fprintf('A è non singolare\n'); 
else
    fprintf('A è singolare\n'); 
end

%% Quadrata
[n,m] = size(A);
if (n == m)
    fprintf('A è una matrice quadrata\n'); 
else
    fprintf('A non è una matrice quadrata\n'); 
end

%% Simmetrica e definita positiva
if isequal(A,A')
	if eig(A) > 0
		fprintf('A è simmetrica definita positiva\n');
	else
		fprintf('A è simmetrica ma non definita positiva\n');
	end
else
	fprintf('A non è simmetrica\n');
end

%% Dominanza diagonale per righe
Adiag = diag(abs(A));
Aout_diag = sum(abs(A),2) - diag(abs(A));
if(Adiag >= Aout_diag)
	if(Adiag > Aout_diag)
		fprintf('La matrice è a dominanza diagonale stretta per righe\n');
	else
		fprintf('La matrice è a dominanza diagonale non stretta per righe\n');
	end
else
	fprintf('La matrice non è a dominanza diagonale per righe\n');
end

%% Dominanza diagonale per colonne
Adiag = diag(abs(A));
Aout_diag = sum(abs(A),1) - diag(abs(A));
if(Adiag >= Aout_diag)
	if(Adiag > Aout_diag)
		fprintf('La matrice è a dominanza diagonale stretta per colonne\n');
	else
		fprintf('La matrice è a dominanza diagonale non stretta per colonne\n');
	end
else
	fprintf('La matrice non è a dominanza diagonale per colonne\n');
end