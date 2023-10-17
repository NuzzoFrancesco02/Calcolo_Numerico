clc
clear

%% caso a forzante nulla
% parametri
K = 100;
L = 20;
N = 20;
n = N-1;

% costruisco la matrice
extradiag = ones (n-1, 1);
maindiag = -2 * ones (n, 1);
A = diag (maindiag, 0) + diag (extradiag, 1) + diag(extradiag, -1);
A = K*A;

% costruisco il termine noto 
t_noto = zeros (n,1);
t_noto(end) = t_noto (end) - K*L;

%risolvo con \ e con thomas
x_mat = A \ t_noto;
[LL, UU, x_thom] = thomas (A, t_noto);

%confronto le soluzioni
figure
plot (abs (x_mat - x_thom), '-o');
title ('differenza fra le soluzioni', 'FontSize', 16)

figure
plot (x_thom, '-o')
title ('posizione', 'FontSize', 16)

figure
plot (x_thom(2:end) - x_thom(1:end-1), '-o')
title ('allungamenti', 'FontSize', 16)

%% caso a forzante non nulla - usiamo rand per generarne una a caso

% costruisco il termine noto 
t_noto=16*rand(n,1);
t_noto(end) = t_noto (end) - K*L;

%risolvo con \ e con thomas
x_mat = A \ t_noto;
[LL, UU, x_thom] = thomas (A, t_noto);

%confronto le soluzioni
figure
plot (abs (x_mat - x_thom), '-o');
title ('differenza fra le soluzioni', 'FontSize', 16)

figure
plot (x_thom, '-o')
title ('posizione', 'FontSize', 16)

figure
plot (x_thom(2:end) - x_thom(1:end-1), '-o')
title ('allungamenti', 'FontSize', 16)

figure;
% plottiamo le forzanti togliendo il termine -K*L dall'ultimo elemento 
% poiche' ha un valore assoluto troppo grande e rende poco significativo
% il grafico
plot ([t_noto(1:end-1); t_noto(end) + K*L], '-o')
title ('forzanti', 'FontSize', 16)
