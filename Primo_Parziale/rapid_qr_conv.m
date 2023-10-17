function mini = rapid_qr_conv(a,b)
l1 = 
l2 = 
l3 =
l4 = 
i = 1;
l = a+0.001:0.001:b+0.01
for j = l
    v(1) = l2/l1;
    v(2) = l3/l2;
    v(3) = l4/l3;
    sol(i) = max(v);
    i = i +1;
end
val = min(sol);
mini = find(sol==val);
l(mini(1))
l(mini(end))