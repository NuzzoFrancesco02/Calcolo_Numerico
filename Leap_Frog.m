function [t,u,v] = Leap_Frog(fun,t_vett,h,y0,w0)
if nargin(fun) ~= 3
    error('Esprimere la funzione come f(t,u,v)')
end
if length(t_vett)~=2
    error('Vettore tempo con più di due elementi');
end
N = (t_vett(2)-t_vett(1))/h;
t = t_vett(1):h:t_vett(2);
u(1) = y0;
v(1) = w0;
for n = 1 : N 
    u(n+1) = u(n) + h*v(n) + ((h^2)/2)*fun(t(n),u(n),v(n));
    phi= @(sol) (v(n) + h*(fun(t(n),u(n),v(n))+fun(t(n+1),u(n+1),sol))/2);
    succ = ptofis(v(n),phi,1e5,1e-5);
    v(n+1) = succ(end);
end
