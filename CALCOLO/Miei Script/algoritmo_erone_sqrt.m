clear
clc

n = 69;
nmax = 1000;
tol = 1e-1;

err = tol + 1;
k = 0;
rk = n;
while(err > tol) && (k < nmax)
	k = k + 1;
	rk(k+1) = 0.5*(rk(k) + n/rk(k));
	err = abs(rk(k+1) - rk(k));
end

figure(1);
clf;
plot(rk);

rk(end)
sqrt(n)