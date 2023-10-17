function [u_n,u] = HN_sym(f,n,y0)
syms t y h s;

N = 0:n;
t_h = N.*h;
u = y0 ;

for i = 1 : n
    s = subs(f,[t y],[t_h(i) u(i)]);
    u = [u simplify(u(i)+h*s)];
    s2 = subs(f,[t y],[t_h(i+1) u(i+1)]);
    u = [u simplify((u(i)+h/2*(s+s2)))];

end
u_n = u(end);

end
