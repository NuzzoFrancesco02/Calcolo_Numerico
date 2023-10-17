function[xk,k] = biseznewton(a,b,fun,dfun,nmax_b,nmax_n,tol)

%Applica il metodo di bisezione per avvicinarsi allo zero e dopo passa al
%metodo di newton

xk = [];
k = [];

%	Bisezione
tol_b = (b-a)/(2^(nmax_b + 1));
[xkb, kb] = bisez(a,b,fun,tol_b);
k = kb;
xk = xkb;

%	Newton
x = xk(end);
[xkn, kn] = newton(x,fun,dfun,nmax_n,tol);
k = k + kn;
xk = [xk; xkn(2:end)];
end