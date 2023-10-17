function fn = twobody(t,y)
[n,m] = size(y);
fn = zeros(n,m);
fn(1) = y(2);
fn(2) = -4*pi^2*y(1)/(y(1)^2+y(3)^2)^(3/2);
fn(3) = y(4);
fn(4) = -4*pi^2*y(3)/(y(1)^2+y(3)^2)^(3/2);
