function [u_n,u] = EA_sym(f,n,y0)
syms t y h s;
N = 0:n;
t_h = N.*h;
u = y0 ;

for i = 1 : n
    s = subs(f,[t y],[t_h(i) u(i)]);
    u = [u simplify(u(i)+h*s)];
end
u_n = u(end);
