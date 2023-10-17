function u_n = EI_sym(f,n,y0)
syms t y h s;

N = 0:n;
t_h = N.*h;
u = y0 ;
for i = 1 : n
    s = subs(f,t,t_h(i+1));
    u = [u simplify(solve(y-u(i)-h*s))];
end
u_n = u(end)