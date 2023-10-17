clc
clear

% parametri
K = 100;
l = 20;

% vettore dei tempi di calcolo
iter = 10;
T = zeros (5, iter);

for i = 1:iter
    % stampo il numero di iterazione  
    disp (strcat ('iter ', num2str(i)))

    %aumento la dimensione
    n = 200*i;

    %costruisco la matrice
    extradiag = ones(n-1, 1);
    maindiag = -2*ones(n, 1);
    A = diag(maindiag, 0) + diag (extradiag, 1) + diag (extradiag, -1);
    A = K*A;
    %costruisco il termine noto
    t_noto = zeros (n, 1);
    %aggiungo la forzante esterna
    t_noto(end) = t_noto(end) - K*l;


    %risolvo con \, LU, thomas, Choleski e memorizzo nella matrice T
    t = tic;
    x_back = A \ t_noto;
    T (1, i) = toc (t);

    t = tic;
    x_back_2 = (-A)\(-t_noto);
    T (2, i) = toc (t);

    t = tic;
    [L, U, P] = lu (A);
    x_lu = U \ (L \ (P*t_noto));
    T (3, i) = toc (t);

    t = tic;
    [L, U, x_thom] = thomas (A, t_noto);
    T (4, i) = toc (t);

    A_chol = -A;
    t_chol = -t_noto;
    t = tic;
    H = chol (A_chol);
    x_chol = H \ (H' \ t_chol);
    T (5, i) = toc (t);
end

figure;
hold on;
dim = 200 * [1:iter];
plot(dim, T(1,:),'-ob', 'linewidth', 2);
plot(dim, T(2,:),'-or', 'linewidth', 2);
plot(dim, T(3,:),'-ok', 'linewidth', 2);
plot(dim, T(4,:),'-om', 'linewidth', 2);
plot(dim, T(5,:),'-og', 'linewidth', 2);

legend('backslash-non-sdp','backslash-sdp', 'lu','thomas','chol', 'Location', 'NorthWest')