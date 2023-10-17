clc
clear
close all

% LABORATORIO 6

%% Punto 1
n    = 50;
A    = diag(4*ones(n,1))     ...
     + diag(-ones(n-1,1),1)  ...
     + diag(-ones(n-2,1),2)  ...
     + diag(-ones(n-1,1),-1) ...
     + diag(-ones(n-2,1),-2);
x0   = zeros(n,1);
b    = 0.2*ones(n,1);
tol  = 1e-5;
nmax = 1e4;

%% Punto 2
disp('Matrice A:')
if (A == A')
   v = eig(A);
   if (v > 0)
      disp('Matrice simmetrica definita positiva')
   else
      disp('Matrice simmetrica ma non definita positiva')
   end
else 
   disp('Matrice non simmetrica')
end

disp('Numero di condizionamento della matrice:')
max(v)/min(v)


%% Punto 4
P = eye(n); % Precondizionatore

disp('Richardson: P = I, alpha = 0.20')
alpha     = 0.2;
B_alpha   = eye(n) - alpha * A; % inv(P)=I
disp('Raggio spettrale:')
max(abs(eig(B_alpha)))
[x02, k02] = richardson(A, b, P, x0, tol, nmax, alpha);


disp('Richardson: P = I, alpha = 0.33')
alpha     = 0.33;
B_alpha   = eye(n) - alpha * A; % inv(P)=I
disp('Raggio spettrale:')
max(abs(eig(B_alpha)))
[x033, k033] = richardson(A, b, P, x0, tol, nmax, alpha);


disp('Richardson: P = I, alpha = OPT')
alpha     = 2/(max(v)+min(v));
B_alpha   = eye(n) - alpha * A; % inv(P)=I
disp('Raggio spettrale:')
max(abs(eig(B_alpha)))
[xopt, kopt] = richardson(A, b, P, x0, tol, nmax, alpha);

% per una formattazione del testo piu' sofisticata proviamo ad utilizzare
% fprintf al posto di disp.

% x_ex = A\b;
% fprintf('Errore: %e\n',norm(x_ex-xopt)/norm(x_ex))
% 
fprintf('\nRichardson non precond. dinamico\n')
[xdin, kdin] = richardson(A, b, P, x0, tol, nmax);



%% Punto 5
P = tril(A); % Precondizionatore

fprintf('\n Richardson: P = tril(A), alpha = 1.00\n')
alpha        = 1.0;
B_alpha      = eye(n) - alpha * inv(P) * A;
fprintf('Raggio spettrale: %f\n', max(abs(eig(B_alpha))))
[x_ri, k_ri] = richardson(A, b, P, x0, tol, nmax, alpha);
[x_gs, k_gs] = gs(A, b, x0, tol, nmax);
fprintf('Gauss-Seidel converge in %d iterazioni \n', k_gs)
fprintf('Scarto tra le soluzioni: %e\n', max(abs(x_ri-x_gs)))


 
%% Punto 6
P = diag(2*ones(n,1))     ...
  + diag(-ones(n-1,1),1)  ...
  + diag(-ones(n-1,1),-1);

fprintf('\nPrecondizionatore P:\n')
if (P == P')
   v = eig(P);
   if (v > 0)
      disp('Matrice simmetrica definita positiva')
   else
      disp('Matrice simmetrica ma non definita positiva')
   end
else 
   disp('Matrice non simmetrica')
end

fprintf('\nRichardson: P = TRIDIAG, alpha = OPT\n')
v         = eig(inv(P)*A);
alpha     = 2/(max(v)+min(v));
B_alpha   = eye(n) - alpha * inv(P) * A;
fprintf('Raggio spettrale: %f\n', max(abs(eig(B_alpha))))
fprintf('Numero di condizionamento: %f\n', cond(inv(P)*A))
[x_tri, k_tri] = richardson(A, b, P, x0, tol, nmax, alpha);


fprintf('\nrichardson: P = TRIDIAG, dinamico\n')
[x, k_tridin] = richardson(A, b, P, x0, tol, nmax);
