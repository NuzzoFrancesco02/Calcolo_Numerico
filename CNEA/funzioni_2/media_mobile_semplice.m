function xt = media_mobile_semplice(x,T,t);
sum = 0;
n = (T-1)/2;
for j = -n:n
    sum = sum + x(t+j);
end
xt = sum/(2*n+1);