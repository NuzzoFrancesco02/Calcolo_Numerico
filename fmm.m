%% x'' = -gamma * x' - omega^2 * x
% y1 = x --> y1' = x' = y2 = f1(t,x,x')=f1(t,y1,y2)
% y2 = x''--> = -gamma*x'-omega^2*x = - gamma*y2 - oemga^2*y1
% = f2(t,y1,y2)
% fn = [f1,f2]
function fn = fmm(t,y)
% y = [y1 y2];
[n,m] = size(y);
fn = zeros(n,m);
gamma = 0.1;
omega2 = 1;
fn(1) = y(2);
fn(2) = -gamma*y(2)-omega2*y(1);
