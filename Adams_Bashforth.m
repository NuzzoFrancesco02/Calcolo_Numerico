function [t,u] = Adams_Bashforth(fun,t_vett,h,y0,u1,u2)
u(1) = y0;
u(2) = u1;
u(3) = u2;
t = t_vett(1):h:t_vett(2);
N = (t_vett(2)-t_vett(1))/h;
for n = 3: N
    u(n+1) = u(n) + h*(23/12*fun(t(n),u(n))-4/3*fun(t(n-1),u(n-1))+5/12*(fun(t(n-2),u(n-2))));
end