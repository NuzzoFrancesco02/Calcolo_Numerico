function [Xk,B,k] = broyden(F,B0,x0,toll,nmax)

% Metodo di Broyden
% Approssima lo zero x della funzione F e la sua jacobiana B

k = 0;
x = x0;
B = B0;
Xk = x;
err = toll + 1;
while(err>toll && k<nmax)
	k = k + 1;
	delta = B \ (-F(x));
	x = x + delta;
	Xk = [Xk,x];
	err = norm(delta);
	B = B + 1/dot(delta,delta)*(F(x) - F(x - delta) - B*delta)*(delta');
end
end