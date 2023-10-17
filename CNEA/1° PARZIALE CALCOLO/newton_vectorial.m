%% [xvect, it, J] = newton_vectorial(F,J,x0,nmax,toll)
% xvect parte da x0 in xvect(:,1)!!!
function [xvect, it, J] = newton_vectorial(F,J,x0,nmax,toll)
if nargin == 4
    toll = 0;
end
err = toll + 1;
x = x0;
it = 0;
xvect = [x0];
while err > toll && it < nmax
    d = J(x) \ (-F(x));
    x = x + d;
    xvect = [xvect x];
    it = it + 1;
    err = norm(F(x));
end
