% n :numero di tempi da iterare
function u_n = CK_sym(f,n,y0)


syms t y h s s2;

N = 0:n;
t_h = N.*h;
u = y0 ;
for i = 1 : n
    s = subs(f,t,t_h(i+1));
    s2 = subs(f,[t y],[t_h(i) u(i)]);

    u = [u simplify(solve(y-u(i)-h/2*(s+s2)))];
   
end
u_n = u(end);
