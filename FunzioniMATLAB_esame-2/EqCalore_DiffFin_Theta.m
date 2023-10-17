function [u, x, t] = EqCalore_DiffFin_Theta(mu, f, a, b, us, ud, g0, ...
	                                        T, h, delta_t, theta)
% [u, x, t] = EqCalore_DiffFin_Theta(f, a, b, us, ud, g0,
%	T, h, delta_t, theta)
% Soluzione del seguente problema per l'equazione del calore:
%   du/dt - mu d2u/dx2 = f       x in (a, b), t in (0, T]
%   u(0,x) = g0(x)            x in (a, b)
%   u(t,a) = us(t)            t in (0, T]
%   u(t,b) = ud(t)            t in (0, T]
%
% INPUT:
%   - mu: coefficiente
%   - f(x,t): forzante
%   - a, b: estremi del dominio
%   - us(t), ud(t): dati di Dirichlet rispettivamente in x = a e x = b
%   - g0(x): dato iniziale
%   - T: tempo finale
%   - h: passo di discretizzazione spaziale
%   - delta_t: passo di discretizzazione temporale
%   - theta: parametro del theta-metodo
% OUTPUT:
%   - u: soluzione numerica; il primo indice determina il nodo
%        nello spazio, il secondo indice l'istante temporale
%        (u(i,j) = soluzione approssimata al nodo x_i, al tempo t_j)
%   - x: vettore dei nodi di discretizzazione
%   - t: vettore dei tempi


% Vettore dei nodi x_j.
Nx = (b - a) / h - 1;
x = linspace(a, b, Nx + 2);

% Vettore dei tempi t^(k).
t = 0:delta_t:T;
Nt = length(t);

% Matrice A.
A = mu / h^2 * (sparse(1:Nx, 1:Nx, 2, Nx, Nx) ...
  + sparse(1:Nx-1, 2:Nx, -1, Nx, Nx) ...
  + sparse(2:Nx, 1:Nx-1, -1, Nx, Nx));

% Matrice A_theta.
A_theta = speye(Nx) + delta_t * theta * A;

% Soluzione approssimata.
u = zeros(length(x), length(t));
u(:,1) = g0(x);

for k = 1:Nt-1
    Fk = f(x(2:end-1), t(k))';
    Fk(1) = Fk(1) + mu * us(t(k)) / h^2;
    Fk(end) = Fk(end) + mu * ud(t(k)) / h^2;

    Fkp1 = f(x(2:end-1), t(k+1))';
    Fkp1(1) = Fkp1(1) + mu * us(t(k+1)) / h^2;
    Fkp1(end) = Fkp1(end) + mu * ud(t(k+1)) / h^2;

    bkp1 = (speye(Nx) - delta_t * (1 - theta) * A) * u(2:end-1,k) ...
         + delta_t * theta * Fkp1 ...
         + delta_t * (1 - theta) * Fk;

    u(2:end-1,k+1) = A_theta \ bkp1;
    u(1,k+1) = us(t(k+1));
    u(end,k+1) = ud(t(k+1));
end

end