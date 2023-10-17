function x = x_float(B,t,m,e,s)
x = 0;
for i = 1 : t
    x = x + m(i)*B^(-i);
end
x = x*((-1)^s)*B^e;