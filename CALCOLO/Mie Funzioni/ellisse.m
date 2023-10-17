function[e, c] = ellisse(x0,y0,a,b);

%Plotta un ellisse (x - x0)/a^2 + (y - y0)/b^2 = 1
% e ne determina eccentricità e, semidistanza folcale c e fuochi

if nargin == 3
	b = a;
end

th = 0:pi/50:2*pi;
xunit = a * cos(th) + x0;
yunit = b * sin(th) + y0;
plot(xunit, yunit);

c = sqrt(abs(a^2 - b^2));

if a >= b
	e = c/a;
	xf1 = x0 + c;
	xf2 = x0 - c;
	fprintf("\nEllisse di eccentricità e = %d, fuochi in F1 = (%d, %d) e F2 = (%d, %d)\n", e, xf1, y0, xf2, y0);
else
	e = c/b;
	yf1 = y0 + c;
	yf2 = y0 - c;
	fprintf("\nEllisse di eccentricità e = %d, fuochi in F1 = (%d, %d) e F2 = (%d, %d)\n", e, x0, yf1, x0, yf2);
end
end