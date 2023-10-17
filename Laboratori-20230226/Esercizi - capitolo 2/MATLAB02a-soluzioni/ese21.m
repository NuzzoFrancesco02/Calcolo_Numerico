clc
clear

% definisco la dimensione
n=1000;

% costruisco le matrici e i termini noti
A = hilb (n);
x_ex = ones (n, 1);
b = A * x_ex;

B = rand (n);
y_ex = ones (n, 1);
c = B * y_ex;

% risolvo i sistemi e calcolo le norme degli errori
x = A \ b;
disp ('---- Matrice A ----')
fprintf ('cond (A) = %5.4e\n', cond(A))
fprintf ('errore relativo = %5.4e\n', (norm (x - x_ex) / norm (x_ex)))
fprintf ('residuo normalizzato = %5.4e\n', (norm (A*x - b) / norm (b)))

y = B \ c;
disp ('---- Matrice B ----')
fprintf ('cond (B) = %5.4e\n', cond(B))
fprintf ('errore relativo = %5.4e\n', (norm (y - y_ex) / norm (y_ex)))
fprintf ('residuo normalizzato = %5.4e\n', (norm (B*y - c) / norm (c)))