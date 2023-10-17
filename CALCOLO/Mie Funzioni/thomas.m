function [L, U, x] = thomas (A, b)
% [L, U, x] = thomas (A, b)
% 
% Calcola la fattorizzazione di Thomas della matrice A e risolve il sistema
% Ax = b
% 
% Parametri in ingresso:
%   - A: matrice tridiagonale quadrata
%   - b: vettore dei termini noti
% 
% Parametri in uscita:
%   - L: la matrice triangolare inferiore della fattorizzazione LU di A
%   - U: la matrice triangolare superiore della fattorizzazione LU di A
%   - x: la soluzione del sistema Ax = b 

%%%%%%%%%%%%%%%%%
% PREPROCESSING %
%%%%%%%%%%%%%%%%%
% dimensione del sistema
N = size (A, 1);
n = length(b);
if size(A,1) ~= n || size(A,2) ~= n || length(x0) ~= n
    error('Le dimensioni non sono compatibili')
end

% crea vettore dei parametri alpha della fattorizzazione di Thomas
alpha = zeros (N, 1);

% crea vettore dei parametri delta della fattorizzazione di Thomas
delta = zeros (N-1, 1);

% estrae sopradiagonale, sottodiagonale e diagonale principale di A
c = diag (A, 1);
e = diag (A, -1);
a = diag (A, 0);


%%%%%%%%%%%%%%%%%%%
% FATTORIZZAZIONE %
%%%%%%%%%%%%%%%%%%%
alpha(1) = a(1);

for i = 2:N
      delta(i-1) = e(i-1) / alpha(i-1);
      alpha(i) = a(i) - delta(i-1) * c(i-1);
end

L = diag (ones (N,1), 0) + diag (delta, -1);
U = diag (alpha, 0) + diag (c, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLUZIONE DEL SISTEMA %
%%%%%%%%%%%%%%%%%%%%%%%%%
y = zeros(N,1);
y(1) = b(1);

for i = 2:N
      y(i) = b(i) - delta(i-1) * y(i-1);
end

x = zeros(N,1);
x(N) = y(N) / alpha(N);

for i = N-1:-1:1
      x(i) = (y(i) - c(i) * x(i+1)) / alpha(i);
end
end