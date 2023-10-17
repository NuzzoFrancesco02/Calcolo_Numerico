function succ = ptofis_sys(x0, phi, nmax, toll)
if size(x0,2)~=1
    x0 = x0';
end

err   = 1 + toll;
it    = 0;
succ  = x0;
xv    = x0;
while (it < nmax && err > toll)
   xn    = phi(xv);
   err   = norm(xn - xv);
   succ = [succ xn];
   it    = it + 1;
   xv    = xn;
end