%% f in funzione di  x e t: f(x,t)
% ATTENZIONE: le colonne partono da t0!!!
% [u, x, t] = EqCalore_DiffFin_Theta(mu, f, a, b, u_s, u_d, g_0, T, h, delta_t,theta)
%
% Delta_t < h^2/(2*mu)
function [u, x, t] = EqCalore_DiffFin_Theta(mu, f, a, b, u_s, u_d, g_0, T, h, delta_t,theta)
if nargin(f) == 1
    error('scrivere f in funzione di  x e t: f(x,t)')
end
Nx = (b-a)/h-1;
x = linspace(a,b,Nx+2);
t = 0:delta_t:T;
Nt = length(t);
A = mu/h^2 * (sparse(1:Nx,1:Nx,2,Nx,Nx) + ...
    sparse(1:Nx-1,2:Nx,-1,Nx,Nx) + ...
    sparse(2:Nx,1:Nx-1,-1,Nx,Nx));
A_theta = speye(Nx)+delta_t*theta*A;

u = zeros(length(x),length(t));
u(:,1) = g_0(x);
for k = 1 : Nt -1
    Fk = f(x(2:end-1),t(k))';
    Fk(1) = Fk(1) + mu/h^2 * u_s(t(k));
    Fk(end) = Fk(end) + mu/h^2 * u_d(t(k));
    
    Fkp1 = f(x(2:end-1),t(k+1))';
    Fkp1(1) = Fkp1(1) + mu/h^2 * u_s(t(k+1));
    Fkp1(end) = Fkp1(end) + mu/h^2* u_d(t(k+1));

    bkp1 = (speye(Nx)-delta_t*(1-theta)*A)*u(2:end-1,k)+...
        delta_t*theta*Fkp1 + ...
        delta_t*(1-theta)*Fk;
    u(2:end-1,k+1) = A_theta\bkp1;
    u(1,k+1) = u_s(t(k+1));
    u(end,k+1) = u_d(t(k+1));
end
