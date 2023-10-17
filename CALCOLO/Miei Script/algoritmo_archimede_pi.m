clear
clc

nmax = 1000;
tol = 1e-6;

an = [sqrt(2)];
bn = [2];

k = 1;
pn = 2^k * an;
qn = 2^k * bn;
err = tol + 1;

while(err > tol) && (k < nmax)
	k = k + 1;
	a = sqrt(2)*sqrt(1 - sqrt(1 - 0.25*(an(k-1)^2)));
	an = [an; a];
	b = an(k)/sqrt(1 - 0.25*(an(k - 1)^2));
	bn = [bn; b];
	pn = [pn; 2^k * an(k)];
	qn = [qn; 2^k * bn(k)];
	err = abs(pn(k) - qn(k));
end

figure(1);
clf;
hold on
plot(an, '--c');
plot(bn, '--b');
plot(pn, 'k');
plot(qn, 'r');

pi
p = pn(end)
q = qn(end)