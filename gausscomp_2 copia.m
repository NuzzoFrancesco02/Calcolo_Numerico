function I = gausscomp_2(a, b, N, f)

h = (b-a)/N;
x = [a:h:b];
xm = [a+h/2:h:b];

x1 = xm- (h/2)*sqrt(15)/5;
x2 = xm;
x3 = xm+ (h/2)*sqrt(15)/5;
y1 = f(x1);
y2 = f(x2);
y3 = f(x3);
I = h/2*(5/9)*sum(y1);
I = I + h/2*(8/9)*sum(y2);
I = I + h/2*(5/9)*sum(y3);
